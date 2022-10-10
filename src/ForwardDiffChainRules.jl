#
# Copyright (c) 2022 The contributors
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

# Adapted from NonconvexUtils.jl with the following license.
#=
MIT License

Copyright (c) 2021 Mohamed Tarek <mohamed82008@gmail.com> and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

module ForwardDiffChainRules

using ChainRulesCore
import ForwardDiff
import MacroTools
using NamedTupleTools

using Requires

include("flatten.jl")

function __init__()
    @require SparseArrays="2f01184e-e22b-5df5-ae63-d93ebab69eaf" begin
        using .SparseArrays: SparseMatrixCSCUnion, SparseVectorUnion, getnzval, getrowval, getcolptr
        include("flatten_SparseArrays.jl")
    end

    @require JuMP="4076af6c-e467-56ae-b986-b466b2749572" begin
        import .JuMP
        include("flatten_JuMP.jl")
    end
end

macro ForwardDiff_frule(sig)
    _fd_frule(sig)
end
export @ForwardDiff_frule

function _fd_frule(sig)
    MacroTools.@capture(sig, f_(x__))
    return quote
        function $(esc(f))($(esc.(x)...))
            f = $(esc(f))
            x = ($(esc.(x)...),)
            flatx, unflattenx = ForwardDiffChainRules.flatten(x)
            CS = length(ForwardDiff.partials(first(flatx)))
            flat_xprimals = ForwardDiff.value.(flatx)
            flat_xpartials = reduce(vcat, transpose.(ForwardDiff.partials.(flatx)))

            xprimals = unflattenx(flat_xprimals)
            xpartials1 = unflattenx(flat_xpartials[:,1])

            yprimals, ypartials1 = ChainRulesCore.frule(
                (NoTangent(), xpartials1...), f, xprimals...,
            )
            flat_yprimals, unflatteny = ForwardDiffChainRules.flatten(yprimals)
            flat_ypartials1, _ = ForwardDiffChainRules.flatten(ypartials1)
            flat_ypartials = hcat(reshape(flat_ypartials1, :, 1), ntuple(Val(CS - 1)) do i
                xpartialsi = unflattenx(flat_xpartials[:, i+1])
                _, ypartialsi = ChainRulesCore.frule((NoTangent(), xpartialsi...), f, xprimals...)
                return ForwardDiffChainRules.flatten(ypartialsi)[1]
            end...)

            T = ForwardDiff.tagtype(eltype(flatx))
            flaty = ForwardDiff.Dual{T}.(
                flat_yprimals, ForwardDiff.Partials.(NTuple{CS}.(eachrow(flat_ypartials))),
            )
            return unflatteny(flaty)
        end
    end
end

function flatten(::Nothing)
    return Float64[], _ -> nothing
end
function flatten(::NoTangent)
    return Float64[], _ -> NoTangent()
end
function flatten(::ZeroTangent)
    return Float64[], _ -> ZeroTangent()
end
function flatten(::Tuple{})
    return Float64[], _ -> ()
end

end # module ForwardDiffChainRules

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
using DifferentiableFlatten

"""
    @ForwardDiff_frule signature mutating

`mutating` indicates whether or not the function is mutating the input argument.

# Example
To define a rule for `LinearAlgebra.exp!`, which is a mutating funciton, we call the macro like this
```julia
@ForwardDiff_frule LinearAlgebra.exp!(A::AbstractMatrix{<:ForwardDiff.Dual}) true
```
"""
macro ForwardDiff_frule(sig, mutating=false)
    _fd_frule(sig; mutating)
end
export @ForwardDiff_frule

using ChainRulesCore
struct ForwardDiffRuleConfig <: RuleConfig{Union{NoReverseMode,HasForwardsMode}} end

function ChainRulesCore.frule_via_ad(
    config::ForwardDiffRuleConfig,
    tangent,
    f,
    f_arg;
    kwargs...,
)
  direcct = frule(tangent, f, f_arg; kwargs...)
  direcct === nothing || return direcct
  T = eltype(f_arg)
  g = t -> f(f_arg + t * tangent[2])
  return g(zero(T)), ForwardDiff.derivative(g, zero(T))
end

const cfg = ForwardDiffRuleConfig()

function _fd_frule(sig; mutating=false)
    if MacroTools.@capture(sig, f_(x__; k__))
        nothing
    else
        MacroTools.@capture(sig, f_(x__))
        k = ()
    end
    return quote
        function $(esc(f))($(esc.(x)...); $(esc.(k)...))
            f = $(esc(f))
            xs = ($(esc.(x)...),)
            ks = (; $(esc.(k)...),)
            flatx, unflattenx = ForwardDiffChainRules.DifferentiableFlatten.flatten(xs)
            CS = length(ForwardDiff.partials(first(flatx)))
            flat_xprimals = ForwardDiff.value.(flatx)
            flat_xpartials = reduce(vcat, transpose.(ForwardDiff.partials.(flatx)))

            xprimals = unflattenx(flat_xprimals)
            xprimals_copy = $mutating ? copy.(xprimals) : xprimals
            xpartials1 = unflattenx(flat_xpartials[:,1])
            yprimals, ypartials1 = ChainRulesCore.frule(
                cfg, (NoTangent(), xpartials1...), f, xprimals_copy...; ks...,
            )
            flat_yprimals, unflatteny = ForwardDiffChainRules.DifferentiableFlatten.flatten(yprimals)
            flat_ypartials1, _ = ForwardDiffChainRules.DifferentiableFlatten.flatten(ypartials1)
            flat_ypartials = hcat(reshape(flat_ypartials1, :, 1), ntuple(Val(CS - 1)) do i
                if $mutating
                    for (xpc, xp) in zip(xprimals_copy, xprimals)
                        xpc .= xp # Update copy
                    end
                end
                xpartialsi = unflattenx(flat_xpartials[:, i+1])
                _, ypartialsi = ChainRulesCore.frule(cfg, (NoTangent(), xpartialsi...), f, xprimals_copy...; ks...)
                return ForwardDiffChainRules.DifferentiableFlatten.flatten(ypartialsi)[1]
            end...)

            T = ForwardDiff.tagtype(eltype(flatx))
            flaty = ForwardDiff.Dual{T}.(
                flat_yprimals, ForwardDiff.Partials.(NTuple{CS}.(eachrow(flat_ypartials))),
            )
            return unflatteny(flaty)
        end
    end
end

end # module ForwardDiffChainRules

#
# Copyright (c) 2022 The contributors
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

using ForwardDiffChainRules

using Test
using ChainRules, ChainRulesCore, ForwardDiff
using LinearAlgebra
import NamedTupleTools: ntfromstruct, structfromnt

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

# TODO tests:
## Multiple outputs
## Struct output
## Functor f - fix first
@testset "ForwardDiff frule" begin
    @eval begin
        f1(x, y) = (x + 2y).^2
        global frule_count = 0
        function ChainRulesCore.frule((_, Δx1, Δx2), ::typeof(f1), x1, x2)
            global frule_count += 1
            println("frule was called")
            return f1(x1, x2), Δx1 + Δx2
        end
        @ForwardDiff_frule f1(x1::ForwardDiff.Dual, x2::ForwardDiff.Dual)
        @ForwardDiff_frule f1(x1::AbstractVector{<:ForwardDiff.Dual}, x2::AbstractVector{<:ForwardDiff.Dual})
        @ForwardDiff_frule f1(x1::AbstractMatrix{<:ForwardDiff.Dual}, x2::AbstractMatrix{<:ForwardDiff.Dual})
        @ForwardDiff_frule LinearAlgebra.exp!(A::AbstractMatrix{<:ForwardDiff.Dual}) true

        f2(x::NamedTuple, y::NamedTuple) = (a = x.a + y.a, b = x.b + y.b)
        f2(x::AbstractVector, y::AbstractVector) = f2.(x, y)
        function ChainRulesCore.frule((_, Δx1, Δx2), ::typeof(f2), x1::NamedTuple, x2::NamedTuple)
            global frule_count += 1
            println("frule was called")
            return f2(x1, x2), (a = Δx1.a + Δx2.a, b = Δx1.b + Δx2.b)
        end
        @ForwardDiff_frule f2(x1::NamedTuple{<:Any, <:Tuple{Vararg{<:ForwardDiff.Dual}}}, x2::NamedTuple{<:Any, <:Tuple{Vararg{<:ForwardDiff.Dual}}})

        struct MyStruct{T, T1, T2}
            a::T1
            b::T2
        end
        MyStruct(a, b) = MyStruct{typeof(a), typeof(a), typeof(b)}(a, b)

        # The @constructor macro takes the type (first) and constructor function (second)
        # The constructor function takes input the fields generated from ntfromstruct (as multiple positional arguments)
        # The ntfromstruct function can be overloaded for your type
        ForwardDiffChainRules.DifferentiableFlatten.@constructor MyStruct MyStruct

        f2(x::MyStruct, y::MyStruct) = MyStruct(x.a + y.a, x.b + y.b)
        function ChainRulesCore.frule((_, Δx1, Δx2), ::typeof(f2), x1::MyStruct, x2::MyStruct)
            global frule_count += 1
            println("frule was called")
            return f2(x1, x2), MyStruct(Δx1.a + Δx2.a, Δx1.b + Δx2.b)
        end
        @ForwardDiff_frule f2(x1::MyStruct{<:ForwardDiff.Dual}, x2::MyStruct{<:ForwardDiff.Dual})
        Base.sum(s::MyStruct) = s.a + s.b

        # I recommend creating your own type to avoid piracy
        _eigvals!(x) = eigvals!(x)
        function ChainRulesCore.frule((_, Δx), ::typeof(_eigvals!), x::Symmetric{<:Real})
            global frule_count += 1
            println("frule was called")
            return frule((NoTangent(), Δx), eigvals!, x)
        end

        # I recommend creating your own type to avoid piracy
        ForwardDiffChainRules.DifferentiableFlatten.@constructor Symmetric Symmetric
        ntfromstruct(a::Symmetric) = (data = a.data,)
        structfromnt(::Type{Symmetric}, x::NamedTuple) = Symmetric(x.data, :U)
        @ForwardDiff_frule _eigvals!(A::Symmetric{<:ForwardDiff.Dual})
    end
    global frule_count = 0
    @testset "2 real inputs - 1 real output" begin
        _f = x -> f1(x[1], x[2])
        _f(rand(2))
        g1 = ForwardDiff.gradient(_f, rand(2))
        @test frule_count == 2
        cfg = ForwardDiff.GradientConfig(_f, rand(2), ForwardDiff.Chunk{2}())
        g2 = ForwardDiff.gradient(_f, rand(2), cfg)
        @test frule_count == 4
        @test g1 == g2
    end
    frule_count = 0
    @testset "2 vector inputs - 1 real output" begin
        _f = x -> sum(f1(x[1:2], x[3:4]))
        g1 = ForwardDiff.gradient(_f, rand(4))
        @test frule_count == 4
        cfg = ForwardDiff.GradientConfig(_f, rand(4), ForwardDiff.Chunk{2}())
        g2 = ForwardDiff.gradient(_f, rand(4), cfg)
        @test frule_count == 8
        @test g1 == g2
    end
    frule_count = 0
    @testset "2 vector inputs - 1 vector output" begin
        _f = x -> f1(x[1:2], x[3:4])
        j1 = ForwardDiff.jacobian(_f, rand(4))
        @test frule_count == 4
        cfg = ForwardDiff.JacobianConfig(_f, rand(4), ForwardDiff.Chunk{2}())
        j2 = ForwardDiff.jacobian(_f, rand(4), cfg)
        @test frule_count == 8
        @test j1 == j2
    end
    frule_count = 0
    @testset "2 matrix inputs - 1 real output" begin
        _f = x -> sum(f1(x[1:2,1:2], x[3:4,3:4]))
        g1 = ForwardDiff.gradient(_f, rand(4, 4))
        @test frule_count == 16
        cfg = ForwardDiff.GradientConfig(_f, rand(4, 4), ForwardDiff.Chunk{2}())
        g2 = ForwardDiff.gradient(_f, rand(4, 4), cfg)
        @test frule_count == 32
        @test g1 == g2
    end
    frule_count = 0
    @testset "2 NamedTuple inputs - 1 real output" begin
        _f = x -> sum(f2((a = x[1], b = x[2]), (a = x[3], b = x[4])))
        g1 = ForwardDiff.gradient(_f, rand(4))
        @test frule_count == 4
        cfg = ForwardDiff.GradientConfig(_f, rand(4), ForwardDiff.Chunk{2}())
        g2 = ForwardDiff.gradient(_f, rand(4))
        @test frule_count == 8
        @test g1 == g2
    end
    frule_count = 0
    @testset "2 struct inputs - 1 real output" begin
        _f = x -> sum(f2(MyStruct(x[1], x[2]), MyStruct(x[3], x[4])))
        g1 = ForwardDiff.gradient(_f, rand(4))
        @test frule_count == 4
        cfg = ForwardDiff.GradientConfig(_f, rand(4), ForwardDiff.Chunk{2}())
        g2 = ForwardDiff.gradient(_f, rand(4))
        @test frule_count == 8
        @test g1 == g2
    end
    frule_count = 0
    @testset "eigvals" begin
        # Gradient of trace
        g = ForwardDiff.gradient(x -> sum(_eigvals!(Symmetric(x))), rand(4, 4))
        @test frule_count == 16
        @test norm(g - I) < 1e-6
    end
    @testset "exp!" begin
        fexp = x -> sum(LinearAlgebra.exp!(copy(x)))
        X = rand(4, 4)
        g = ForwardDiff.gradient(fexp, X)
        g2 = FiniteDifferences.grad(central_fdm(5, 1), fexp, X)[1]
        @test norm(g-g2) < 1e-4
    end
    @testset "kwargs" begin
        fkwarg(x1, x2; a = 2.0) = x1 * x2 * a
        frule_count = 0
        function ChainRulesCore.frule((_, Δx1, Δx2), ::typeof(fkwarg), x1::Real, x2::Real; a)
            global frule_count += 1
            println("frule was called")
            return fkwarg(x1, x2; a), a * x2 * Δx1 + a * x1 * Δx2
        end
        @ForwardDiff_frule fkwarg(x1::ForwardDiff.Dual, x2::ForwardDiff.Dual; kwargs...)
        @test ForwardDiff.gradient(x -> fkwarg(x[1], x[2], a = 3.0), [1.0, 2.0]) == [6, 3]
        @test frule_count == 2
    end
end

# ForwardDiffChainRules.jl
[![Run Tests](https://github.com/ThummeTo/ForwardDiffChainRules.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/ThummeTo/ForwardDiffChainRules.jl/actions/workflows/Test.yml)
[![Coverage](https://codecov.io/gh/ThummeTo/ForwardDiffChainRules.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ThummeTo/ForwardDiffChainRules.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

## What is ForwardDiffChainRules.jl?
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) does not support `frule`s from [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) by default. As a result, if you are creating custom AD-rules and want support for the most common AD-tools in Julia, you need to define three differentiation rules: 
- `ChainRulesCore.rrule`
- `ChainRulesCore.frule`
- an additional dispatch of your function for values of type `ForwardDiff.Dual`
Technically, the last two candidates aim both for forward sensitivities, so include the same differentiation rules. This is redundant code and an error prone coding task... and not necessary anymore! 

[ForwardDiffChainRules.jl](https://github.com/ThummeTo/ForwardDiffChainRules.jl) allows you to re-use the differentiation code defined in an existing `ChainRulesCore.frule` with only a few lines of code and without re-coding your differentiation rules.

## How can I use ForwardDiffChainRules.jl?
1\. Open a Julia-REPL, switch to package mode using `]`, activate your preferred environment.

2\. Install [*ForwardDiffChainRules.jl*](https://github.com/ThummeTo/ForwardDiffChainRules.jl):
```julia-repl
(@v1.6) pkg> add ForwardDiffChainRules
```

3\. If you want to check that everything works correctly, you can run the tests bundled with [*ForwardDiffChainRules.jl*](https://github.com/ThummeTo/ForwardDiffChainRules.jl):
```julia-repl
(@v1.6) pkg> test ForwardDiffChainRules
```

4\. Have a look inside the [examples folder](https://github.com/ThummeTo/ForwardDiffChainRules.jl/tree/main/examples).

## How can I add a dispatch for ForwardDiff based on an existing `frule`? 
```julia
using ForwardDiffChainRules

function f1(x1, x2)
   # do whatever you want to do in your function
   return (x + 2y).^2
end

# define your frule for function f1 as usual
function ChainRulesCore.frule((_, Δx1, Δx2), ::typeof(f1), x1, x2)
   # this could be any code you want of course
   return f1(x1, x2), Δx1 + Δx2
end

# create a ForwardDiff-dispatch for scalar type `x1` and `x2`
@ForwardDiff_frule f1(x1::ForwardDiff.Dual, x2::ForwardDiff.Dual)

# create a ForwardDiff-dispatch for vector type `x1` and `x2`
@ForwardDiff_frule f1(x1::AbstractVector{<:ForwardDiff.Dual}, x2::AbstractVector{<:ForwardDiff.Dual})

# create a ForwardDiff-dispatch for matrix type `x1` and `x2`
@ForwardDiff_frule f1(x1::AbstractMatrix{<:ForwardDiff.Dual}, x2::AbstractMatrix{<:ForwardDiff.Dual})
```

## Acknowledgement
This package is based on code from Mohamed Tarek (@mohamed82008) in his package [NonconvexUtils.jl](https://github.com/JuliaNonconvex/NonconvexUtils.jl). The initial discussion started on [discourse.julialang.org](https://discourse.julialang.org/t/chainrulescore-and-forwarddiff/61705). With the aim of providing this functionality as light-weigth as possible, this package was created.

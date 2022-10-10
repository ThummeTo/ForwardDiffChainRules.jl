# ForwardDiffChainRules.jl
[![Run Tests](https://github.com/ThummeTo/ForwardDiffChainRules.jl/actions/workflows/Test.yml/badge.svg)](https://github.com/ThummeTo/ForwardDiffChainRules.jl/actions/workflows/Test.yml)
[![Coverage](https://codecov.io/gh/ThummeTo/ForwardDiffChainRules.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ThummeTo/ForwardDiffChainRules.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

## What is ForwardDiffChainRules.jl?
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) does not support frules from [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) by default. As a result, if you are creating custom AD-rules and want support for the most common AD-tools in Julia, you need to define three differentiation rules: `ChainRulesCore.rrule`, `ChainRulesCore.frule` and a dispatch of your function for values of type `ForwardDiff.Dual`. Technically, the last two candidates aim both for forward sensitivities, so include the same differentiation operations. This is redundant code and error prone... and not necessary anymore! 

[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) allows you to re-use the differentiation code defined in an existing `ChainRulesCore.frule` with only a few lines of code and without re-coding your differnetiation operations.

## How can I use ForwardDiffChainRules.jl?
1\. Open a Julia-REPL, switch to package mode using `]`, activate your preferred environment.

2\. Install [*ForwardDiffChainRules.jl.jl*](https://github.com/ThummeTo/ForwardDiffChainRules.jl.jl):
```julia-repl
(@v1.6) pkg> add ForwardDiffChainRules.jl
```

3\. If you want to check that everything works correctly, you can run the tests bundled with [*ForwardDiffChainRules.jl.jl*](https://github.com/ThummeTo/ForwardDiffChainRules.jl.jl):
```julia-repl
(@v1.6) pkg> test ForwardDiffChainRules.jl
```

4\. Have a look inside the [examples folder](https://github.com/ThummeTo/ForwardDiffChainRules.jl.jl/tree/main/examples).

## How can I add a dispatch for ForwardDiff based on an existing `frule`? 
```julia
using ForwardDiffChainRules
```

## Acknowledgement
This package is based on code from Mohamed Tarek (@mohamed82008) in his package [NonconvexUtils.jl](https://github.com/JuliaNonconvex/NonconvexUtils.jl). The initial discussion started on [discourse.julialang.org](https://discourse.julialang.org/t/chainrulescore-and-forwarddiff/61705). With the aim of providing this functionality in a package as ligth-weigth as possible, this package was created.
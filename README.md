# GinzburgLandau

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://maximilian-gelbrecht.github.io/GinzburgLandau.jl/dev/)
[![Build Status](https://github.com/maximilian-gelbrecht/GinzburgLandau.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/maximilian-gelbrecht/GinzburgLandau.jl/actions/workflows/CI.yml?query=branch%3Amain)

Two-Dimensional discretization of the Complex Ginzburg Landau Equation $\partial_t u = (1 + i\alpha)\Delta u + u - (1+i\beta)|u|^2u$. The package does not include a solver, but needs to used in conjuction with e.g. `DifferentialEquations.jl`. Includes a finite differences and a pseudo-spectral discretization. Right now, only periodic boundary conditions are implemented. For more, see the [Documentation]((https://maximilian-gelbrecht.github.io/GinzburgLandau.jl/dev/))

## Usage 

The implementation models the phase state as a 1D flattened array that can be reshaped easily into a 2D array on the grid. 

### Finite Differences 

```julia 
using GinzburgLandau, OrdinaryDiffEq, Plots

n = 128 # number of grid points per side
L = 192 # domain size of each side

α = 2.0
β = -1.0

g = GinzburgLandau.GridFD(range(0,L,length=n),range(0,L,length=n))
Δ = GinzburgLandau.Laplacian2DPeriodic(g)
u0 = GinzburgLandau.initial_conditions(g)
p = [α, β]
tspan = (0.,200.)

prob = ODEProblem(GinzburgLandau.cgle_fd!, u0, tspan, [α, β, Δ])

sol = solve(prob, Tsit5())

ts = tspan[1]:0.25:tspan[2]
anim = @animate for t in ts
    heatmap(reshape(abs.(sol(t)), g), clim=(0.,1.5))
end

gif(anim, "cgle2d.gif")
```

![cgle2d.gif](cgle2d.gif)

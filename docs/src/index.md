```@meta
CurrentModule = GinzburgLandau
```

# GinzburgLandau

Two-Dimensional discretization of the Complex Ginzburg Landau Equation $\partial_t u = (1 + i\alpha)\Delta u + u - (1+i\beta)|u|^2u$. The package does not include a solver, but needs to used in conjuction with e.g. `DifferentialEquations.jl`. Includes a finite differences and a pseudo-spectral discretization. Right now, only periodic boundary conditions are implemented. 

## Usage 

### Finite Differences 

The implementation models the phase state as a 1D flattened array that can be reshaped easily into a 2D array on the grid. A discretization of the domain is performed and a 2D finite difference operator constructed. 

```julia 
using GinzburgLandau, OrdinaryDiffEq, Plots

n = 128 # number of grid points per side
L = 192 # domain size of each side

α = 2.0
β = -1.0

g = GinzburgLandau.GridFD(range(0,L,length=n),range(0,L,length=n))
Δ = GinzburgLandau.Laplacian2DPeriodic(g)
u0 = GinzburgLandau.initial_conditions(g)
p = [α, β, Δ]
tspan = (0.,200.)

prob = ODEProblem(GinzburgLandau.cgle_fd!, u0, tspan, p)

sol = solve(prob, Tsit5())

ts = tspan[1]:0.25:tspan[2]
anim = @animate for t in ts
    heatmap(reshape(abs.(sol(t)), g), clim=(0.,1.5))
end

gif(anim, "cgle2d.gif")
```
## Pseudospectral 

Here, the phase state is actually 2D all the time and the equation is solve in the Fourier domain in which the Laplace operator simply evalautes to ``k^2`` with the wavenumber ``k``. However, as the reaction term is nonlinear, it is not easily computable in the Fourier domain, so that a pseudospectral approach, i.e. transforming back and forward, is applied. 

```julia
using OrdinaryDiffEq, GinzburgLandau, StatsBase, Plots,FFTW 

n = 36
L = 75
α = 2.0
β = -1.0

g = GinzburgLandau.GridSP(range(0,L,length=n),range(0,L,length=n))
p = GinzburgLandau.CGLESPeudoSpectralPars(g, α, β)
u0 = p.FT * GinzburgLandau.initial_conditions(g)
tspan = (0.,200.)

prob = ODEProblem(GinzburgLandau.cgle_ps!, u0, tspan, p)
sol = solve(prob, Tsit5())
sol_array = Array(sol(0.:0.5:200.))

# transform back 
sol_real = similar(sol_array)
for i ∈ 1:size(sol_array,3)
    sol_real[:,:,i] = p.iFT * sol_array[:,:,i]
end

# plot
anim = @animate for i ∈ 1:size(sol_real, 3)
    heatmap(abs.(sol_real[:,:,i]), clim=(0.,1.5))
end

gif(anim, "cgle2d-ps.gif")
```
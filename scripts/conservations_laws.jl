using GinzburgLandau, OrdinaryDiffEq, Plots

n = 128 # number of grid points per side
L = 192 # domain size of each side

α = 2.0
β = -1.0

g = GinzburgLandau.GridFD(range(0,L,length=n),range(0,L,length=n))
Δ = GinzburgLandau.Laplacian2DPeriodic(g)
u0 = GinzburgLandau.initial_conditions(g)
tspan = (0.,200.)

prob = ODEProblem(GinzburgLandau.cgle_fd!, u0, tspan, [α, β, Δ])

sol = solve(prob, Tsit5(), reltol=1e-6)

sol_matrix = reshape(Array(sol(50.:1:200.)),n,n,:)

PLOTS = false
if PLOTS 
    anim = @animate for it in axes(sol_matrix,3)
        heatmap(abs.(sol_matrix[:,:,it]), clim=(0.,1.5))
    end
    gif(anim, "cgle2d.gif")
end

function impulse_power(u, h=1) 
    q² = abs2.(u)
    P = (h*h) * reshape(sum(q², dims=(1,2)),:)
    return P 
end 



impulse_power(u, g::AbstractGrid) = impulse_power(u, g.h)
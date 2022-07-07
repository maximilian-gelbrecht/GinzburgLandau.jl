n = 36
L = 75
α = 2.0
β = -1.0

g = GinzburgLandau.GridFD(range(0,L,length=n),range(0,L,length=n))
Δ = GinzburgLandau.Laplacian2DPeriodic(g)
u0 = GinzburgLandau.initial_conditions(g)
p = [α, β]
tspan = (0.,200.)

prob = ODEProblem(GinzburgLandau.cgle_fd!, u0, tspan, [α, β, Δ])

sol = solve(prob, Tsit5())

sol_array = Array(sol(50.:1:200.))

# test that mean and variance stays rougly constant after some transient
@test sum(0.75 .< mean(abs.(sol_array), dims=1)[:] .< 0.85) == size(sol_array,2)
@test sum(0.045 .< var(abs.(sol_array), dims=1)[:] .< 0.075) == size(sol_array,2)

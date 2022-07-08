using FFTW 

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

sol_array = Array(sol(50.:0.5:200.))


# transform back 
sol_real = similar(sol_array)
for i ∈ 1:size(sol_array,3)
    sol_real[:,:,i] = p.iFT * sol_array[:,:,i]
end

# test that mean and variance stays rougly constant after some transient
@test sum(0.7 .< mean(abs.(sol_real), dims=(1,2))[:] .< 0.9) == size(sol_real,3)
@test sum(0.045 .< var(abs.(sol_real), dims=(1,2))[:] .< 0.075) == size(sol_real,3)



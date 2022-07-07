begin
    using FFTW
    using Plots
    using OrdinaryDiffEq
end

begin
    # simple meshgrid to get the coordinate axes both in real and spectral space
    function meshgrid(x::AbstractArray{T,1}) where T<:Number
        N = length(x)
        x1 = zeros(T, (N,N))
        y1 = zeros(T, (N,N))

        for i=1:N
            x1[i,:] = x
        end
        y1 = collect(x1')

        return x1,y1
    end

    struct CGLE2d_Pars{T,S,U,V}
        LinOp::AbstractArray{T,2} # Linear Differential Operator
        β::S
        FT::U # Fourier Transform Operator
        iFT::V # Inverse Fourier Transform Operator
    end

    function CGLE2d_Pars(n::Int, L::Int, α::T, β::T, u₀::AbstractArray{T,2}) where T<:Number
        @assert iseven(n)
        n2 = Int(n/2)

        # the grid in spectral domain
        k = [0:n2-1; 0; -n2+1:-1] .*(2π/L);
        k² = k.*k
        k²[n2 + 1] = (n2*(2π/L))^2
        kx², ky² = meshgrid(k²)
        ksum = kx² .+ ky²

        # the linear operator
        LinFac = 1.0 .- (ksum .* (1.0 .+ α))

        FT = FFTW.plan_fft(u₀)
        iFT = FFTW.plan_ifft(FT*u₀)

        CGLE2d_Pars(LinFac, β, FT, iFT)
    end
end

begin
    n = 36
    L = 75
    α = 2.0im
    β = -1.0im
    #u0 = 0.1*(rand(ComplexF64, (n,n)) .- 0.5)
    dx = n/L
    x,y = meshgrid(range(-L/2, stop=L/2, length=n))
    u₀ = 0.01*(rand(ComplexF64, (n,n)) .- 0.5)
    pars = CGLE2d_Pars(n, L, α, β, u₀)
    Fu₀ = pars.FT*u₀

    non_lin(x, β) = x - (1 + β)*abs(x)^2*x
end

# finite difference version
function cgle_fd!(du, u, p, t)
    α, β = p
    du .= (1 .+ α)*lap*u + non_lin.(u, β)
end

# pseudo-spectral
function cgle!(du, u, p, t)
    invu = p.iFT*u
    du .= p.LinOp.*u .- p.FT*( (1.0 .+p.β) .* (abs.(invu)).^2 .* invu)
end

t_end = 100. # short integration time for initial experiment (the animation is pretty slow to render)
prob = ODEProblem(cgle!, Fu₀, (0.,t_end), pars)
t_start = 50. #throw away transient

println("solving....")
if LOAD_DATA
    @load "ginsburg-data.jld2" dat
else
    @time u = solve(prob, Tsit5())
    #throw away transient
    t = t_start:0.1:t_end
    dat = zeros(eltype(u₀), n, n, length(t))
    for it=1:length(t)
        dat[:,:,it] = pars.iFT*u(t[it])
    end
    if SAVE_DATA
        @save "ginsburg-data.jld2" dat
    end
end

N_t = length(50.:0.1:t_end)

if PLOT
    Plots.pyplot()
    anim = Plots.@animate for i in eachindex(t)
        Plots.heatmap(abs.(dat[:,:,i]), clim=(-1.5,1.5))
    end

    gif(anim, "clge2d-small-25.gif")
end

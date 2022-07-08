using AbstractFFTs 

"""
    GridSP

2D dicretization grid for the pseudo spectral approach

# Initializiation 

    GridFD(x::AbstractRange, y::AbstractRange)

# Fields

* `x` x coordinate, a ``n``-long range
* `y` y coordinate, a ``n``-long range
* `h` grid spacing 
* `n` number of grid points along one side
* `N` total number of grid points 
* `k` wavenumber 
* `k²` squared wavenumber 
* `k²_x_grid` meshgrid for `k_x`, this is a ``n\times n`` matrix 
* `k²_y_grid` meshgrid for `k_y`
"""
struct GridSP <: AbstractGrid
    x # this is a n_x long range
    y
    h # spacing 
    n # number of grid points along one side
    N # total number of grid points 
    k # wave number 
    k²
    k²_x_grid 
    k²_y_grid
end 

function GridSP(x::AbstractRange, y::AbstractRange)
    h_x = abs(x[2] - x[1])
    h_y = abs(y[2] - y[1])
    n_x = length(x)
    n_y = length(y)
    @assert n_x == n_y
    @assert h_x == h_y 

    N = n_x * n_y 

    L = abs(x[end] - x[1])

    @assert iseven(n_x)
    n2 = Int(n_x/2)
    k = [0:n2-1; 0; -n2+1:-1] .*(2π/L);
    k² = k.*k
    k²[n2 + 1] = (n2*(2π/L))^2

    k²_x_grid = ones(n_x) *  k²'
    k²_y_grid = k² * ones(n_y)'
    
    GridSP(x, y, h_x, n_x, N, k, k², k²_x_grid, k²_y_grid)
end

initial_conditions(g::GridSP) = initial_conditions(g.n)

"""
    CGLESPeudoSpectralPars{T,S,U,V}

Parameter struct for Pseudospectral approach

    CGLESPeudoSpectralPars(g::GridSP, α, β) 

# Fields 

* `LinOp::AbstractArray{T,2}` Linear Differential Operator
* `α::S`
* `β::S`
* `FT::U` Fourier Transform Operator
* `iFT::V` Inverse Fourier Transform Operator
"""
struct CGLESPeudoSpectralPars{T,S,U,V}
    LinOp::AbstractArray{T,2} # Linear Differential Operator
    α::S
    β::S
    FT::U # Fourier Transform Operator
    iFT::V # Inverse Fourier Transform Operator
end

function CGLESPeudoSpectralPars(g::GridSP, α, β) 
        
    ksum = g.k²_x_grid .+ g.k²_y_grid

    # the linear operator
    LinOp = 1.0 .- (ksum .* complex(1.0,α))

    u₀ = initial_conditions(g)

    FT = plan_fft(u₀)
    iFT = plan_ifft(FT*u₀)

    CGLESPeudoSpectralPars(LinOp, α,  β, FT, iFT)
end

function cgle_ps!(du, u, p, t)
    LinOp, α, β, FT, iFT = p.LinOp, p.α, p.β, p.FT, p.iFT
    invu = iFT*u
    du .= LinOp.*u .- FT*(complex(1.0,p.β) .* (abs.(invu)).^2 .* invu)
end 


using SparseArrays, LinearAlgebra

"""
    GridFD

2D dicretization grid for the finite difference approach

# Initializiation 

    GridFD(x::AbstractRange, y::AbstractRange)

# Fields

* `x` x coordinate, a ``n``-long range
* `y` y coordinate, a ``n``-long range
* `h` grid spacing 
* `n` number of grid points along one side
* `N` total number of grid points 
* `x_grid` meshgrid for `x`, this is a ``n\times n`` matrix 
* `y_grid` meshgrid for `y`
"""
struct GridFD <: AbstractGrid
    x # this is a n_x long range
    y
    h # spacing 
    n # number of grid points along one side
    N # total number of grid points 
    x_grid # this is a n_x x n_y matrix 
    y_grid
end 

function GridFD(x::AbstractRange, y::AbstractRange)
    h_x = abs(x[2] - x[1])
    h_y = abs(y[2] - y[1])
    n_x = length(x)
    n_y = length(y)
    @assert n_x == n_y
    @assert h_x == h_y 

    N = n_x * n_y 
    x_grid = ones(n_x) * x'
    y_grid = y * ones(n_y)'
    
    GridFD(x, y, h_x, n_x, N, x_grid, y_grid)
end

"""
    Laplacian2DPeriodic(g::GridFD)

Returns 2d finite diffence laplace operator for a domain of side length `n` with spacing `dx` as a sparse matrix. Periodic boundary condtions are applied
"""
function Laplacian2DPeriodic(g::GridFD)
    n = g.n
    N = g.N 
    
    # these are the neighbours in x-direction
    # however this will include points that are not actually neighbours in the grid (n,n+1),(n+1,n)
    # we thus will manually delete those form the array and replace them with the correct neighours 
    # at the left and right boundary in the end 
    M = diagm(0=>-4*ones(N), 1=>ones(N-1), -1=>ones(N-1))
    for i=1:(n-1) 
        M[i*n,i*n + 1] = 0
        M[i*n + 1,i*n] = 0
    end

    # these are the neighbours in y-direction
    M += diagm(n=>ones(N-n), -n=>ones(N-n))

    # these are the neighours of the grid points at the upper boundary
    M += diagm(N-n => ones(n))

    # these are the neighbours of the grid points at the lower boundary 
    M += diagm(-N+n => ones(n))
    
    for i=1:n:N # loop over left and right boundary points
       M[i,i+(n-1)] = 1
       M[i+(n-1),i] = 1
    end 

    sparse(M./(g.h^2))
end

"""
    Laplacian2DPeriodic(n::Integer, L)

Returns 2d finite diffence laplace operator for a domain of side length `L` and `n` grid points as a sparse matrix. Periodic boundary condtions are applied
"""
Laplacian2DPeriodic(n::Integer, L) = Laplacian2DPeriodic(GridFD(range(0, L, length=n), range(0, L, length=n)))

function cgle_fd!(du, u, p, t)
    α, β, Δ = p
    du .= complex(1,α).*(Δ*u) .+ reaction(u, β)
end
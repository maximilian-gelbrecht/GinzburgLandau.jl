# functions that both dicretizations need 
abstract type AbstractGrid end 
import Base.reshape

"""
    initial_conditions(n::Integer)

Returns (random) initial conditions that are commonly used, `0.01*(rand(ComplexF64, (n,n)) .- (0.5 + 0.5im))`
"""
initial_conditions(n::Integer) = 0.01*(rand(ComplexF64, (n,n)) .- (0.5 + 0.5im))
initial_conditions(g::AbstractGrid) = reshape(initial_conditions(g.n), g)

"""
    reshape(x::AbstractArray, g::AbstractGrid)

Reshapes between 2D array and a flattened 1D array.
"""
Base.reshape(x::AbstractArray{T,2}, g::AbstractGrid) where T = reshape(x, :)
Base.reshape(x::AbstractArray{T,1}, g::AbstractGrid) where T = reshape(x, g.n, g.n)

"""
    reaction(x, β)

The reaction term of the equation: ``u - (1+i\beta)|u|^2u``
"""
reaction(x, β) = x .- complex(1,β) .* abs.(x).^2 .* x


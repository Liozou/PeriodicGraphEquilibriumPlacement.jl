# Main functions

using PeriodicGraphs, Graphs, StaticArrays

"""
    equilibrium(g::PeriodicGraph)

Return an equilibrium placement for the vertices of the graph, defined as a list
of positions such that each vertex is at the barycentre of its neighbors.

The returned equilibrium placement is such that the first vertex of the graph
is at the origin of the space.
"""
function equilibrium(g::PeriodicGraph{N}) where N
    n = nv(g)
    iszero(n) && return Matrix{Rational{Int64}}(undef, N, 0)
    Y = Matrix{Int}(undef, n, N)
    A = spzeros(Int, n, n)
    neigh = Vector{Int}(undef, n)
    offset = SizedVector{N,Int}(undef)
    for i in 1:n
        neigh .= 0
        offset .= 0
        count = 0
        for k in neighbors(g, i)
            k.v == i && continue
            count += 1
            neigh[k.v] += 1
            offset .-= k.ofs
        end
        Y[i,:] .= offset
        A[i,:] .= neigh
        A[i,i] = -count
    end

    return dixon_solve(Val(N), A[2:end,2:end], Y[2:end,:])
end

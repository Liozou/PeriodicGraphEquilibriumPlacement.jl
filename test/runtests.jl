using PeriodicGraphEquilibriumPlacement
using Graphs, PeriodicGraphs, SparseArrays
using Test, Random

function is_at_mean_position(g::PeriodicGraph{D}, poss, i) where D
    neighs = neighbors(g, i)
    (@view poss[:,i]) == sum((@view poss[:,v]) .+ o for (v, o) in neighs) .// length(neighs)
end

function all_at_mean_position(g, poss=equilibrium(g))
    all(x -> is_at_mean_position(g, poss, x), 1:nv(g))
end

@testset "equilibrium" begin
    pcu = PeriodicGraph("3 1 1 0 0 1 1 1 0 1 0 1 1 1 0 0");
    @test equilibrium(pcu) == reshape(Rational{Int64}[0//1; 0//1; 0//1], 3, 1)
    @test all_at_mean_position(pcu)

    dia = PeriodicGraph("3 1 2 0 0 0 1 2 0 0 1 1 2 0 1 0 1 2 1 0 0");
    @test all_at_mean_position(dia)

    srs = PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 2 3 1 0 0 2 4 0 -1 0 3 4 0 0 -1");
    @test all_at_mean_position(srs)

    afy = PeriodicGraph("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 1 5 0 0 0 2 6 0 0 0 2 7 0 0 0 2 8 0 0 0 3 7 0 0 0 3 9 0 0 0 3 10 0 0 0 4 8 0 0 0 4 9 0 0 0 4 11 0 0 0 5 7 0 1 0 5 8 0 1 0 5 9 0 1 0 6 12 0 0 0 6 13 0 0 0 6 14 0 0 0 7 14 0 0 0 8 15 0 0 0 9 12 -1 0 1 10 12 -1 0 1 10 13 -1 0 1 10 15 0 0 1 11 13 -1 0 0 11 14 -1 0 0 11 15 0 0 0 12 16 0 0 0 13 16 0 1 0 14 16 0 0 0 15 16 -1 0 0")
    @test all_at_mean_position(afy)
end

@testset "dixon_solve and rational_solve" begin
    for N in 1:5
        for n in 1:5
            A = sparse(Int.(rand(Int8, n, n)))
            Y = Int.(rand(Int8, n, N))
            result = dixon_solve(Val(N), A, Y)
            @test A*result == Y
            @test result == rational_solve(Val(N), A, Y)
        end
    end
end

@testset "dixon_solve edge cases" begin
    A = Int.(rand(Int8, 7, 7))
    Y = Int.(rand(Int8, 7, 3))
    for (i,p) in ((5, 2147483647), (7, 2147483629), (4, 2147483587))
        v = zeros(Int, 7)
        v[i] = p
        A[:,i] = A[i,:] = v
        @test A*(dixon_solve(Val(3), sparse(A), Y)) == Y
    end
end

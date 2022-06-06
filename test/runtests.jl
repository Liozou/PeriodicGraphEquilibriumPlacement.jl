using PeriodicGraphEquilibriumPlacement, PeriodicGraphs, Graphs
using Test

function is_at_mean_position(g::PeriodicGraph{D}, poss, i) where D
    neighs = neighbors(g, i)
    (@view poss[:,i]) == sum((@view poss[:,v]) .+ o for (v, o) in neighs) .// length(neighs)
end

function all_at_mean_position(g, poss=equilibrium(g))
    all(x -> is_at_mean_position(g, poss, x), 1:nv(g))
end

@testset "PeriodicGraphEquilibriumPlacement.jl" begin
    pcu = PeriodicGraph("3 1 1 0 0 1 1 1 0 1 0 1 1 1 0 0");
    @test equilibrium(pcu) == reshape(Rational{Int64}[0//1; 0//1; 0//1], 3, 1)
    @test all_at_mean_position(pcu)

    dia = PeriodicGraph("3 1 2 0 0 0 1 2 0 0 1 1 2 0 1 0 1 2 1 0 0");
    @test all_at_mean_position(dia)
end

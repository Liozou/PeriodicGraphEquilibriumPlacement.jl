using PeriodicGraphEquilibriumPlacement
using Graphs, PeriodicGraphs, SparseArrays
using Test, Random, LinearAlgebra

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

    wrong = PeriodicGraph2D(PeriodicGraph("1  1 2 0  1 3 0  5 2 0  5 3 0  2 3 1"))
    @test_throws ErrorException equilibrium(wrong)
end

@testset "dixon_solve and rational_solve" begin
    for N in 1:5
        for n in 1:5
            _A = Int.(rand(Int8, n, n))
            A = sparse(Int.(_A) .+ 300*LinearAlgebra.I(n)) # ensure the system is invertible
            Y = Int.(rand(Int8, n, N))
            result = dixon_solve(Val(N), A, Y)
            @test A*result == Y
            @test result == rational_solve(Val(N), A, Y)
        end
    end
    A = sparse([208 72 887 687 946 263 905 943 131 183 256 606 613 854 914; 582 279 104 638 1 272 130 214 910 813 376 376 202 362 384; 786 526 26 700 84 417 430 253 894 41 895 207 620 918 163; 753 98 421 556 839 665 861 678 614 245 548 186 831 774 642; 834 257 952 786 485 529 66 833 619 258 886 901 488 55 87; 917 536 333 316 295 528 645 777 236 4 247 641 411 101 262; 394 87 654 584 354 858 361 570 990 326 279 348 984 623 251; 774 33 215 61 577 492 634 769 755 577 176 989 964 379 104; 176 21 886 253 198 886 545 136 958 175 311 646 954 730 927; 654 216 757 545 975 391 56 480 258 861 639 126 53 295 739; 29 619 956 175 693 567 734 415 645 704 818 291 678 557 973; 964 112 555 991 366 297 105 482 26 848 221 687 50 53 119; 155 938 591 60 807 704 157 929 289 128 846 971 31 658 530; 736 711 941 890 816 473 171 60 299 167 434 819 582 913 856; 668 14 20 336 336 823 651 164 468 737 181 828 192 93 58])
    Y = [15*j + i for i in 1:15, j in -1:1]
    result = dixon_solve(Val(3), A, Y)
    @test A*result == Y
    @test result == rational_solve(Val(3), A, Y)
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

    A2 = sparse([1, 2, 4, 1, 2, 4, 1, 2, 4], [1, 1, 1, 2, 2, 2, 4, 4, 4], [-3, 1, 1, 1, -3, 1, 1, 1, -2], 4, 4)
    Y2 = [1 1; 0 2; 1 -1; 0 0];
    @test isempty(rational_solve(Val(2), A2, Y2))
    @test isempty(dixon_solve(Val(2), A2, Y2))
end

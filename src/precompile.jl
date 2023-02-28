using PeriodicGraphEquilibriumPlacement, PeriodicGraphs, LinearAlgebra, SparseArrays

const primes = (2147483647, 2147483629, 2147483587)
const modulo{P} = PeriodicGraphEquilibriumPlacement.Modulo{P, Int32}

macro enforce(expr) # strong @assert
    msg = string(expr)
    return :($(esc(expr)) ? $(nothing) : throw(AssertionError($msg)))
end


function _precompile_()
    # SparseArrays
    @enforce precompile(Tuple{typeof(*),SparseArrays.SparseMatrixCSC{Int, Int},Matrix{Rational{BigInt}}})
    @enforce precompile(Tuple{typeof(Base.copyto_unaliased!),IndexCartesian,SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},IndexLinear,Vector{Int}})
    @enforce precompile(Tuple{typeof(Base.mightalias),SubArray{Int, 1, SparseArrays.SparseMatrixCSC{Int, Int}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, false},Vector{Int}})
    @enforce precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{BigInt},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{BigInt},Bool,Bool})
    @enforce precompile(Tuple{typeof(SparseArrays.dimlub),Vector{Int}})
    @enforce precompile(Tuple{typeof(SparseArrays.findnz),SparseArrays.SparseMatrixCSC{Int, Int}})
    @enforce precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{Int},Int,Type})
    @enforce precompile(Tuple{typeof(SparseArrays.spzeros),Type{Int},Type{Int},Int,Int})
    @enforce precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{Int, Int},UnitRange{Int},UnitRange{Int}})

    # Modulos.jl
    for P in primes
        @enforce precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{modulo{P}},Int,Int,Function})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{modulo{P}},Int,Int,Function})
        @enforce precompile(Tuple{typeof(-),modulo{P},modulo{P}})
        @enforce precompile(Tuple{typeof(/),Int,modulo{P}})
        @enforce precompile(Tuple{typeof(==),Matrix{modulo{P}},Matrix{Int}})
        @enforce precompile(Tuple{typeof(Base._unsafe_copyto!),Matrix{modulo{P}},Int,Matrix{Int},Int,Int})
        @enforce precompile(Tuple{typeof(Base._unsafe_getindex),IndexLinear,Matrix{modulo{P}},Int,Base.Slice{Base.OneTo{Int}}})
        @enforce precompile(Tuple{typeof(Base.copyto_unaliased!),IndexLinear,SubArray{modulo{P}, 1, Matrix{modulo{P}}, Tuple{Int, Base.Slice{Base.OneTo{Int}}}, true},IndexLinear,Vector{modulo{P}}})
        @enforce precompile(Tuple{typeof(LinearAlgebra.mul!),Matrix{modulo{P}},SparseArrays.SparseMatrixCSC{Int, Int},Matrix{modulo{P}},Bool,Bool})
        @enforce precompile(Tuple{typeof(SparseArrays._setindex_scalar!),SparseArrays.SparseMatrixCSC{modulo{P}, Int},Int,Int,Int})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse!),Vector{Int},Vector{Int},Vector{modulo{P}},Int,Int,typeof(+),Vector{Int},Vector{Int},Vector{Int},Vector{modulo{P}},Vector{Int},Vector{Int},Vector{modulo{P}}})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse),Vector{Int},Vector{Int},Vector{modulo{P}},Int,Int,Function})
        @enforce precompile(Tuple{typeof(SparseArrays.sparse_check_length),String,Vector{modulo{P}},Int,Type})
        @enforce precompile(Tuple{typeof(getindex),SparseArrays.SparseMatrixCSC{modulo{P}, Int},UnitRange{Int},UnitRange{Int}})
    end

    # solver.jl
    @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.rational_lu!), SparseMatrixCSC{Rational{BigInt},Int}, Vector{Int}, Bool})
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.rational_lu!), SparseMatrixCSC{Rational{BigInt},Int}, Vector{Int}})
    for P in primes
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.rational_lu!), SparseMatrixCSC{modulo{P},Int}, Vector{Int}, Bool})
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.rational_lu!), SparseMatrixCSC{modulo{P},Int}, Vector{Int}})
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.rational_lu), SparseMatrixCSC{modulo{P},Int}, Bool, Type{modulo{P}}})
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.rational_lu), SparseMatrixCSC{modulo{P},Int}, Bool, Type{Rational{BigInt}}})
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.rational_lu), SparseMatrixCSC{Int,Int}, Bool})
    end
    for Ti in (Rational{BigInt}, (modulo{P} for P in primes)...)
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.forward_substitution!), SparseMatrixCSC{Ti,Int}, Matrix{Ti}})
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.backward_substitution!), SparseMatrixCSC{Ti,Int}, Matrix{Ti}})
    end
    @static if VERSION < v"1.8-"
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.linsolve!), LU{Rational{BigInt},SparseMatrixCSC{Rational{BigInt},Int}}, Matrix{Rational{BigInt}}})
        for P in primes
            @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.linsolve!), LU{modulo{P},SparseMatrixCSC{modulo{P},Int}}, Matrix{Int}})
        end
    else
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.linsolve!), LU{Rational{BigInt},SparseMatrixCSC{Rational{BigInt},Int},Base.OneTo{Int}}, Matrix{Rational{BigInt}}})
        for P in primes
            @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.linsolve!), LU{modulo{P},SparseMatrixCSC{modulo{P},Int},Base.OneTo{Int}}, Matrix{Int}})
        end
    end
    for T in (Int64, Int128, BigInt)
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.copyuntil), Int, Matrix{Rational{T}}, Type{Rational{T}}})
        @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement._inner_dixon_p!), Vector{Int}, Matrix{Rational{T}}, BigInt, Matrix{BigInt}, BigInt, BigInt})
    end
    for N in 1:3
        @enforce precompile(Tuple{typeof(rational_solve), Val{N}, SparseMatrixCSC{Int,Int}, Matrix{Int}})
        for P in primes
            @static if VERSION < v"1.8-"
                @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.dixon_p), Val{N}, SparseMatrixCSC{Int,Int}, LU{modulo{P},SparseMatrixCSC{modulo{P},Int}}, Matrix{Int}})
            else
                @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement.dixon_p), Val{N}, SparseMatrixCSC{Int,Int}, LU{modulo{P},SparseMatrixCSC{modulo{P},Int},Base.OneTo{Int}}, Matrix{Int}})
            end
        end
        @enforce precompile(Tuple{typeof(dixon_solve), Val{N}, SparseMatrixCSC{Int,Int}, Matrix{Int}})
    end

    # embeddings.jl
    for N in 1:3
        for T in (Rational{Int64}, Rational{Int128}, Rational{BigInt})
            @enforce precompile(Tuple{typeof(PeriodicGraphEquilibriumPlacement._catzeros), Val{N}, Adjoint{T,Matrix{T}}})
        end
        @enforce precompile(Tuple{typeof(equilibrium), PeriodicGraph{N}})
    end
end

_precompile_()

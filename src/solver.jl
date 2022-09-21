## Specialized solver for sparse symmetric integer linear systems with exact rational solution

using .Modulos
import Base.GMP: MPZ

using Base: OneTo
using LinearAlgebra: BlasInt, checknonsingular, LU, tril!, triu!, issuccess, norm
using SparseArrays
using SparseArrays: getcolptr

using BigRationals


function compat_lu(::Type{Tf}, B, maxvec, info) where Tf
    @static if VERSION < v"1.8-"
        LU{Tf,SparseMatrixCSC{Tf,Int}}(dropzeros!(B), Vector{BlasInt}(1:maxvec), convert(BlasInt, info))
    else
        LU{Tf,SparseMatrixCSC{Tf,Int},OneTo{Int}}(Tf.(B), 1:maxvec, info)
    end
end

function compat_lu_convert(::Type{Tf}, B, maxvec, info) where Tf
    @static if VERSION < v"1.8-"
        LU{Tf,SparseMatrixCSC{Tf,Int}}(Tf.(dropzeros!(B)), Vector{BlasInt}(1:maxvec), convert(BlasInt, info))
    else
        LU{Tf,SparseMatrixCSC{Tf,Int},OneTo{Int}}(Tf.(B), 1:maxvec, info)
    end
end

function rational_lu!(B::SparseMatrixCSC, col_offset::Vector{Int}, check=true)
    Tf = eltype(B)
    m, n = size(B)
    minmn = min(m, n)
    info = 0
    @inbounds begin
        for k in 1:minmn
            ipiv = getcolptr(B)[k] + col_offset[k]
            piv = nonzeros(B)[ipiv]
            if iszero(piv)
                check && checknonsingular(k-1, Val(false)) # TODO update with Pivot
                return compat_lu(Tf, B, minmn, k-1)
            end
            Bkkinv = inv(piv)
            nzB = nonzeros(B)
            @simd ivdep for i in ipiv+1:getcolptr(B)[k+1]-1
                nzB[i] *= Bkkinv
            end
            for j in k+1:n
                r1 = getcolptr(B)[j]
                r2 = getcolptr(B)[j+1]-1
                r = searchsortedfirst(rowvals(B), k, r1, r2, Base.Forward)
                ((r > r2) || (rowvals(B)[r] != k)) && continue
                Bkj = nonzeros(B)[r]
                for i in ipiv+1:getcolptr(B)[k+1]-1
                    Bik = nonzeros(B)[i]
                    l = i - ipiv

                    while rowvals(B)[l+r] < rowvals(B)[i]
                        r += 1
                    end
                    nonzeros(B)[l+r] -= Bik * Bkj
                end
            end
        end
    end
    check && checknonsingular(info, Val(false))
    return compat_lu(Tf, B, minmn, info)
end

# function lu!(B::SparseMatrixCSC{<:Rational}, ::Val{Pivot} = Val(false);
#                            col_offset, check::Bool = true) where Pivot
function rational_lu!(B::SparseMatrixCSC{BigRational}, col_offset::Vector{Int}, check::Bool=true)
    Tf = Rational{BigInt}
    m, n = size(B)
    minmn = min(m, n)
    info = 0
    Bkkinv = BigRational()
    tmp = BigRational()
    @inbounds begin
        for k in 1:minmn
            ipiv = getcolptr(B)[k] + col_offset[k]
            piv = nonzeros(B)[ipiv]
            if iszero(piv)
                check && checknonsingular(k-1, Val(false)) # TODO update with Pivot
                return compat_lu_convert(Tf, B, minmn, k-1)
            end
            BigRationals.MPQ.inv!(Bkkinv, piv)
            @simd for i in ipiv+1:getcolptr(B)[k+1]-1
                BigRationals.MPQ.mul!(nonzeros(B)[i], Bkkinv)
            end
            for j in k+1:n
                r1 = getcolptr(B)[j]
                r2 = getcolptr(B)[j+1]-1
                r = searchsortedfirst(rowvals(B), k, r1, r2, Base.Forward)
                ((r > r2) || (rowvals(B)[r] != k)) && continue
                Bkj = nonzeros(B)[r]
                for i in ipiv+1:getcolptr(B)[k+1]-1
                    Bik = nonzeros(B)[i]
                    l = i - ipiv

                    while rowvals(B)[l+r] < rowvals(B)[i]
                        r += 1
                    end
                    # Base.GMP.MPZ.mul!(tmp, Bik, Bkj)
                    # Base.GMP.MPZ.sub!(nonzeros(B)[l+r], tmp)
                    BigRationals.MPQ.mul!(tmp, Bik, Bkj)
                    BigRationals.MPQ.sub!(nonzeros(B)[l+r], tmp)
                end
            end
        end
    end
    check && checknonsingular(info, Val(false))
    return compat_lu_convert(Tf, B, minmn, info)
end

# function lu(A::SparseMatrixCSC{<:Rational}, pivot::Union{Val{false}, Val{true}} = Val(false); check::Bool = true)
function rational_lu(A::SparseMatrixCSC, check::Bool=true, ::Type{Ti}=BigRational) where {Ti}
    Tf = Ti == BigRational ? Rational{BigInt} : Ti

    Base.require_one_based_indexing(A)
    _I, _J, _V = findnz(A)
    I, J, V = issorted(_J) ? (_I, _J, _V) : begin
        _indices = sortperm(_J)
        @inbounds (_I[_indices], _J[_indices], _V[_indices])
    end
    # @inbounds if !issorted(_J)
    #     indices = sortperm(J)
    #     I = I[indices]; J = J[indices]; V = V[indices]
    # end
    isempty(J) && return compat_lu_convert(Tf, A, 0, 0)
    m, n = size(A)
    minmn = min(m, n)
    if J[1] != 1 || I[1] != 1
        check && checknonsingular(1, Val(false)) # TODO update with Pivot
        # return LU{eltype(A), typeof(A)}(A, collect(1:minmn), convert(BlasInt, 1))
        return compat_lu_convert(Tf, A, minmn, 1)
    end

    col_offset = zeros(Int, minmn) # for each col, index of the pivot element
    idx_cols = [[I[i] for i in getcolptr(A)[col+1]-1:-1:getcolptr(A)[col]] for col in 1:minmn]
    # For each column, indices of the non-zeros elements
    in_idx_colscol = falses(n)
    for col in 2:minmn
        sort!(idx_cols[col-1]; rev=true)
        # All idx_cols[x] are sorted by decreasing order for x < col
        # @show idx_cols[col]
        idx_colscol = idx_cols[col]
        in_idx_colscol[idx_colscol] .= true
        for row_j in idx_colscol
            row_j >= col && continue
            col_offset[col] += 1
            # @show idx_cols[row_j]
            idx_colsj = idx_cols[row_j]
            sizcol = length(idx_colscol)
            for row_i in idx_colsj
                if row_i ≤ row_j
                    break # Because the row_i are sorted in decreasing order
                end
                if !in_idx_colscol[row_i]
                    push!(idx_colscol, row_i)
                    in_idx_colscol[row_i] = true
                end
            end
            countadd = length(idx_colscol) - sizcol
            if countadd > 0
                siz = length(I)
                resize!(I, siz + countadd)
                resize!(J, siz + countadd)
                resize!(V, siz + countadd)
                for i in 1:countadd
                    row_i = idx_colscol[sizcol+i]
                    _idx = siz + i
                    J[_idx] = col
                    I[_idx] = row_i
                    V[_idx] = 0
                end
            end
        end
        in_idx_colscol[idx_colscol] .= false
    end
    B = sparse(I, J, Ti.(V)) # TODO update with Pivot
    rational_lu!(B, col_offset, check)
end


function forward_substitution!(L::SparseMatrixCSC, b)
    _, n = size(L)
    _, m = size(b)
    @inbounds for col in 1:n
        k = getcolptr(L)[col]
        if rowvals(L)[k] != col && col != 0
            return false
        end
        invnzLk = inv(nonzeros(L)[k])
        x = invnzLk .* b[col,:]
        b[col,:] .= x
        for i in (k+1):getcolptr(L)[col+1]-1
            nzLi = nonzeros(L)[i]
            rvLi = rowvals(L)[i]
            @simd ivdep for j in 1:m
                b[rvLi,j] -= nzLi*x[j]
            end
        end
    end
    true
end

function backward_substitution!(U::SparseMatrixCSC, b)
    _, n = size(U)
    _, m = size(b)
    @inbounds for col in n:-1:1
        k = getcolptr(U)[col+1]-1
        if rowvals(U)[k] != col && col != 0
            return false
        end
        invnzUk = inv(nonzeros(U)[k])
        x = invnzUk .* b[col,:]
        b[col,:] .= x
        for i in getcolptr(U)[col]:(k-1)
            nzUi = nonzeros(U)[i]
            rvUi = rowvals(U)[i]
            @simd ivdep for j in 1:m
                b[rvUi,j] -= nzUi*x[j]
            end
        end
    end
    true
end

function linsolve!(F::LU, B::Matrix)
    TFB = typeof(oneunit(eltype(B)) / oneunit(eltype(F)))
    BB = similar(B, TFB, size(B))
    copyto!(BB, B)
    m, n = size(F)
    minmn = min(m,n)
    L = tril!(getfield(F, :factors)[1:m, 1:minmn])
    for i = 1:minmn; L[i,i] = 1; end
    forward_substitution!(L, BB) || return BB, false
    x = triu!(getfield(F, :factors)[1:minmn, 1:n])
    backward_substitution!(x, BB) || return BB, false
    return BB, true
end

"""
    rational_solve(::Val{N}, A::SparseMatrixCSC{Int,Int}, Y::Matrix{Int}) where N

Fallback solver for [`dixon_solve`](@ref) which performs an LU decomposition followed by
forward and backward substitutions.

In general, it is slower than [`dixon_solve`](@ref).
"""
function rational_solve(::Val{N}, A::SparseMatrixCSC{Int,Int}, Y::Matrix{Int}) where N
    B = rational_lu(A, false)
    if !issuccess(B)
        error("Singular exception while equilibrating. Is the graph connected?")
    end
    Z, check = linsolve!(B, Rational{BigInt}.(Y))
    check || error("Singular exception on substitution. Please report this error by opening an issue.")
    return Rational{Int128}.(Z)
    # Rational{Int64} is not enough for tep for instance.
end



function copyuntil(j, oldZ, ::Type{T}) where T
    Z = similar(oldZ, T)
    for i in eachindex(Z)
        i == j && return Z
        Z[i] = oldZ[i]
    end
    error("Invalid failure of _inner_dixon_p!. Please report this error by opening an issue.")
    return Z # Does not matter but just in case for type stability
end


function _inner_dixon_p!(indices::Vector{Int}, Z::Matrix{Rational{T}}, h::BigInt,
                         x̄::Matrix{BigInt}, sqh::BigInt, tmp::BigInt) where T
    while !isempty(indices)
        j = pop!(indices)
        ua = MPZ.set(h)
        ub = deepcopy(@inbounds x̄[j])
        va = Int128(0)
        vb = Int128(1)
        k = 0
        while ub >= sqh
            k += 1
            # cpua = deepcopy(ua)
            # cpub = deepcopy(ub)
            MPZ.tdiv_qr!(tmp, ua, ua, ub)
            ua, ub = ub, ua
            # @assert tmp == cpua ÷ cpub
            # @assert ua == cpub
            # @assert ub == cpua - tmp * cpub
            # cpuc = deepcopy(va)
            if typemin(Clong) < vb < typemax(Clong)
                MPZ.mul_si!(tmp, vb % Clong)
            else
                tmp *= vb
            end
            flag = signbit(va)
            va = abs(va)
            if va < typemax(Culong)
                if flag
                    MPZ.sub_ui!(tmp, va)
                else
                    MPZ.add_ui!(tmp, va)
                end
                va, vb = vb, Int128(tmp)
            else
                va, vb = vb, va + tmp
            end
            #= or replace all of the above since if typemin(Clong) < ... by
            MPZ.mul!(tmp, vb)
            MPZ.add!(tmp, va)
            va, vb, tmp = vb, tmp, va
            =#
            # @assert vb == cpuc + tmp * va
        end

        uv::Tuple{T,T} = if T === BigInt
            Base.divgcd(ub, vb)
        else
            ud, vd = Base.divgcd(ub, vb)
            m = typemin(T)
            M = typemax(T)
            if !(m < ud < M && m < vd < M)
                push!(indices, j)
                return false
            end
            (ud % T, vd % T)
        end

        @inbounds Z[j] = (-2*isodd(k)+1) * Base.checked_den(uv[1], uv[2])
        # @assert mod((-1)^isodd(k) * ub, h) == mod(vb * x̄[j], h)
    end
    return true
end

function dixon_p(::Val{N}, A::SparseMatrixCSC{Int,Int}, C::LU{Modulo{p,Int32}}, Y::Matrix{Int}) where {N,p}
    λs = [norm(x) for x in eachcol(A)]
    append!(λs, norm(x) for x in eachcol(Y))
    partialsort!(λs, N)
    for _ in 1:N
        popfirst!(λs)
    end
    δ::BigFloat = prod(BigFloat, λs; init=one(BigFloat))
    m = ceil(Int, 2*log(δ / (MathConstants.φ - 1))/log(p))
    # @assert m ≥ 1
    B = copy(Y)
    Z::Union{Matrix{Rational{Int64}},Matrix{Rational{Int128}},Matrix{Rational{BigInt}}} = similar(Y, Rational{Int64})
    BB, check = linsolve!(C, B)
    check || return Z, false
    x̄ = BigInt.(BB)
    X = copy(x̄)
    # @assert A * Modulo{p,Int32}.(X) == B
    h = one(BigInt) # = p^i
    tmp = BigInt()
    for _ in 1:m-1
        MPZ.mul_si!(h, p)
        B .= (B .- A*Integer.(X)) .÷ p
        BB2, check2 = linsolve!(C, B)
        check2 || return Z, false
        X .= Integer.(BB2)
        # @assert A * Modulo{p,Int32}.(X) == B
        @inbounds for j in eachindex(x̄)
            MPZ.mul!(tmp, X[j], h)
            MPZ.add!(x̄[j], tmp)
        end
    end
    MPZ.mul_si!(h, p) # h = p^m
    # @assert mod.(A * x̄, h) == mod.(Y, h)
    sqh = MPZ.sqrt(h) # h = p^{m/2}
    indices = collect(reverse(eachindex(Z)))
    success = _inner_dixon_p!(indices, Z, h, x̄, sqh, tmp)
    if !success
        Z = copyuntil(first(indices), Z, Rational{Int128})
        success = _inner_dixon_p!(indices, Z, h, x̄, sqh, tmp)
        if !success
            Z = copyuntil(first(indices), Z, Rational{BigInt})
            success = _inner_dixon_p!(indices, Z, h, x̄, sqh, tmp)
            # @assert success
        end
    end

    # @assert eltype(Y).(A * big.(Z)) == Y
    return Z, true
end

@static if VERSION < v"1.8-"
    const typeofB = Union{
        LU{Modulo{2147483647,Int32},SparseMatrixCSC{Modulo{2147483647,Int32},Int}},
        LU{Modulo{2147483629,Int32},SparseMatrixCSC{Modulo{2147483629,Int32},Int}},
        LU{Modulo{2147483587,Int32},SparseMatrixCSC{Modulo{2147483587,Int32},Int}}
    }
else
    const typeofB = Union{
        LU{Modulo{2147483647,Int32},SparseMatrixCSC{Modulo{2147483647,Int32},Int},OneTo{Int}},
        LU{Modulo{2147483629,Int32},SparseMatrixCSC{Modulo{2147483629,Int32},Int},OneTo{Int}},
        LU{Modulo{2147483587,Int32},SparseMatrixCSC{Modulo{2147483587,Int32},Int},OneTo{Int}}
    }
end

function try_modulo(::Val{N}, A, Y, ::Type{Modulo{p,T}}) where {N,p,T}
    B::typeofB = rational_lu(A, false, Modulo{p,Int32})
    issuccess(B) || Matrix{Rational{Int64}}(undef, 0, 0), false
    return dixon_p(Val(N), A, B, Y)
end

"""
    dixon_solve(::Val{N}, A::SparseMatrixCSC{Int,Int}, Y::Matrix{Int}) where N

Specialized solver for the linear system `A*X = Y` where `A` is a sparse integer `n×n`
matrix and `Y` is a dense integer `n×N` matrix, using Dixon's method.

Return `X` as either a `Matrix{Rational{Int64}}`, a `Matrix{Rational{Int128}}` or a
`Matrix{Rational{BigInt}}`, whichever smallest type can hold all its values.
"""
function dixon_solve(::Val{N}, A::SparseMatrixCSC{Int,Int}, Y::Matrix{Int}) where N
    # @show time_ns()
    Z, success = try_modulo(Val(N), A, Y, Modulo{2147483647,Int32})
    success && @goto ret
    Z, success = try_modulo(Val(N), A, Y, Modulo{2147483629,Int32})
    success && @goto ret
    Z, success = try_modulo(Val(N), A, Y, Modulo{2147483587,Int32})
    success && @goto ret
    # The probability of this being required is *extremely* low
    return rational_solve(Val(N), A, Y)
    @label ret
    return Z
end

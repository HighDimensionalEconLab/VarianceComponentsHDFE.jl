# Computes the LDLinv, along with the adjaceny matrix, both need to be passed
# around to create a solver for a parallelized thread
function computeLDLinv(sddm)
    a,d = adj(sddm)
    a1 = extendMatrix(a,d)
    ldli = approxchol_lap_pc(a1)
    la = lap(a1)
    return ldli, la
end



# Wrap the LDLinv object as a preallocated linear operator
function approxcholOperator(ldli::LDLinv{Tind,Tval},buff::Vector{Tval}) where {Tind,Tval}
    prod = @closure rhs -> LDLsolver!(buff,ldli,rhs)
    return PreallocatedLinearOperator{Tval}(length(ldli.d), length(ldli.d), true, true, prod, nothing, nothing)
end

# Wrap the Algebraic Multigrid object as preallocated linear operator
function AmgOperator(ml::AlgebraicMultigrid.MultiLevel, buff::Vector{Tval}) where {Tval}
    prod =  @closure rhs -> AlgebraicMultigrid.solve!(buff,ml,rhs)
    return PreallocatedLinearOperator{Tval}(length(buff), length(buff), true, true, prod, nothing, nothing)
end

# Compute a solver for a grounded system (SDDM matrix) with a PreallocatedLinearOperator, and an adjacency matrix
function approxcholSolver(P::PreallocatedLinearOperator, la::AbstractArray; tol::Real=1e-6, maxits=300, verbose=0)

    tol_=tol
    maxits_=maxits
    verbose_=verbose

    f = function(b;tol=tol_, maxits=maxits_, verbose=verbose_)
        xaug = Krylov.cg(la,[b; -sum(b)] .- mean([b; -sum(b)]), M=P, rtol = tol, itmax=maxits, verbose=0)[1]
        xaug .= xaug .- xaug[end]
        return xaug[1:end-1]
    end

    return f

end

# Compute a solver for a grounded system (SDDM matrix) with a LDLinv object, and an adjacency matrix
function approxcholSolver(ldli::LDLinv, la::AbstractArray; tol::Real=1e-6, maxits=300, verbose=0)

    buff = zeros(length(ldli.d))
    P = approxcholOperator(ldli,buff)

    tol_=tol
    maxits_=maxits
    verbose_=verbose

    f = function(b;tol=tol_, maxits=maxits_, verbose=verbose_)
        xaug = Krylov.cg(la,[b; -sum(b)] .- mean([b; -sum(b)]), M=P, rtol = tol, itmax=maxits, verbose=0)[1]
        xaug .= xaug .- xaug[end]
        return xaug[1:end-1]
    end

    return f

end

# Compute a solver for a grounded system (SDDM matrix)
function approxcholSolver(sddm::AbstractArray; tol::Real=1e-6, maxits=300, verbose=0)

    # Compute LDLinv object and adjacency matrix
    ldli,la = computeLDLinv(sddm)

    # Wrap the LDLinv object as a preallocated linear operator
    buff = zeros(length(ldli.d))
    P = approxcholOperator(ldli,buff)

    tol_=tol
    maxits_=maxits
    verbose_=verbose

    f = function(b;tol=tol_, maxits=maxits_, verbose=verbose_)
        xaug = Krylov.cg(la,[b; -sum(b)] .- mean([b; -sum(b)]), M=P, rtol = tol, itmax=maxits, verbose=0)[1]
        xaug .= xaug .- xaug[end]
        return xaug[1:end-1]
    end

    return f

end



########## UTILS FUNCTION BROUGHT FROM FIXEDEFFECTSMODELS.JL 
function invsym!(X::Symmetric; setzeros = false, diagonal = 1:size(X, 2))
    tols = max.(diag(X), 1)
    buffer = zeros(size(X, 1))
    for j in diagonal
        d = X[j,j]
        if setzeros && abs(d) < tols[j] * sqrt(eps())
            X.data[1:j,j] .= 0
            X.data[j,(j+1):end] .= 0
        else
            # used to mimic SAS; now similar to SweepOperators
            copy!(buffer, view(X, :, j))
            Symmetric(BLAS.syrk!('U', 'N', -1/d, buffer, one(eltype(X)), X.data))
            rmul!(buffer, 1 / d)
            @views copy!(X.data[1:j-1,j], buffer[1:j-1])        
            @views copy!(X.data[j, j+1:end], buffer[j+1:end])   
            X[j,j] = - 1 / d
        end
        if setzeros && j == 1
            tols = max.(diag(X), 1)
        end
    end
    return X
end


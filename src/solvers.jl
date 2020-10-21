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



# Compute a solver for a grounded system (SDDM matrix) with a PreallocatedLinearOperator, and an adjacency matrix
function approxcholSolver(P::PreallocatedLinearOperator, la::AbstractArray; tol::Real=1e-6, maxits=300, verbose=false)

    tol_=tol
    maxits_=maxits
    verbose_=verbose

    f = function(b;tol=tol_, maxits=maxits_, verbose=verbose_)
        xaug = Krylov.cg(la,[b; -sum(b)] .- mean([b; -sum(b)]), M=P, rtol = tol, itmax=maxits, verbose=verbose)[1]
        xaug .= xaug .- xaug[end]
        return xaug[1:end-1]
    end

    return f

end

# Compute a solver for a grounded system (SDDM matrix) with a LDLinv object, and an adjacency matrix
function approxcholSolver(ldli::LDLinv, la::AbstractArray; tol::Real=1e-6, maxits=300, verbose=false)

    buff = zeros(length(ldli.d))
    P = approxcholOperator(ldli,buff)

    tol_=tol
    maxits_=maxits
    verbose_=verbose

    f = function(b;tol=tol_, maxits=maxits_, verbose=verbose_)
        xaug = Krylov.cg(la,[b; -sum(b)] .- mean([b; -sum(b)]), M=P, rtol = tol, itmax=maxits, verbose=verbose)[1]
        xaug .= xaug .- xaug[end]
        return xaug[1:end-1]
    end

    return f

end

# Compute a solver for a grounded system (SDDM matrix)
function approxcholSolver(sddm::AbstractArray; tol::Real=1e-6, maxits=300, verbose=false)

    # Compute LDLinv object and adjacency matrix
    ldli,la = computeLDLinv(sddm)

    # Wrap the LDLinv object as a preallocated linear operator
    buff = zeros(length(ldli.d))
    P = approxcholOperator(ldli,buff)

    tol_=tol
    maxits_=maxits
    verbose_=verbose

    f = function(b;tol=tol_, maxits=maxits_, verbose=verbose_)
        xaug = Krylov.cg(la,[b; -sum(b)] .- mean([b; -sum(b)]), M=P, rtol = tol, itmax=maxits, verbose=verbose)[1]
        xaug .= xaug .- xaug[end]
        return xaug[1:end-1]
    end

    return f

end

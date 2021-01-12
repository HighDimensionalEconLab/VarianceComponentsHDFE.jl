
## Dependencies

#To install Laplacians you need to put in Julia's Prompt "add Laplacians#master"
#For some reason Pkg.add() doesn't work with this.
#Pkg.add("MATLAB")
#using Laplacians
#using Pkg
#include(string(Pkg.dir("Laplacians") , "/src/matlabSolvers.jl"))
getlagged(x) = [NaN; x[1:(end - 1)]]

using DocStringExtensions


#Defining types and structures
abstract type AbstractLLSAlgorithm end
abstract type AbstractGraphLLSAlgorithm <: AbstractLLSAlgorithm end  # i.e. requires graph laplacian
struct CMGPreconditionedLLS <: AbstractGraphLLSAlgorithm  end
struct AMGPreconditionedLLS <: AbstractGraphLLSAlgorithm  end
struct DirectLLS <: AbstractLLSAlgorithm  end  # no graph required

abstract type AbstractLeverageAlgorithm end

"""
$(TYPEDEF)

Data type to pass to VCHDFESettings type, to indicate Exact algorithm
"""
struct ExactAlgorithm <: AbstractLeverageAlgorithm  end

"""
$(TYPEDEF)

Data type to pass to VCHDFESettings type, to indicate JLA algorithm

### Fields

* `num_simulations`: number of simulations in estimation. If num_simulations = 0, defaults to 100 * log(#total fixed effect)"
"""
@with_kw struct JLAAlgorithm <: AbstractLeverageAlgorithm
    num_simulations::Int64 = 0
end

"""
$(TYPEDEF)

The VCHDFESettings type is to pass information to methods regarding which algorithm to use. 

### Fields

* `cg_maxiter`: maximum number of iterations (default = 300)
* `leverage_algorithm`: which type of algorithm to use (default = ExactAlgorithm())
* `first_id_effects`: includes first id effects. At this version it is required to include the first_id_effects. (default = true)
* `cov_effects`: includes covariance of first-second id effects. At this version it is required to include the cov_effects. (default = true)
* `print_level`: prints the state of the program in std output. If print_level = 0, the app prints nothing in the std output. (default = 1)
* `first_id_display_small`: name of the first id in lower cases (default = person)
* `first_id_display`: name of the first id (default = Person)
* `second_id_display_small`: name of the second id in lower cases (default = firm)
* `second_id_display`: name of the second id (default = Firm)
* `outcome_id_display_small`: name of the observation id in lower cases (default = wage)
* `outcome_id_display`: name of the observation id (default = Wage)

"""
@with_kw struct VCHDFESettings{LeverageAlgorithm}
    cg_maxiter::Int64 = 300
    leverage_algorithm::LeverageAlgorithm = ExactAlgorithm()
    leave_out_level::String = "obs"
    first_id_effects::Bool = true
    cov_effects::Bool = true
    print_level::Int64 = 1
    first_id_display_small::String = "person"
    first_id_display::String = "Person"
    second_id_display_small::String = "firm"
    second_id_display::String = "Firm"
    outcome_id_display_small::String = "wage"
    outcome_id_display::String = "Wage"
end

## FUNCTIONS FROM THE VECTORIZEDROUTINE.JL PACKAGE

function accumarray(subs, val, sz=(maximum(subs),))
    A = zeros(eltype(val), sz...)
    for i = 1:length(val)
        @inbounds A[subs[i]] += val[i]
    end
    A
end

function accumarray(subs, val::Number, sz=(maximum(subs),))
    A = zeros(eltype(val), sz...)
    for i = 1:length(subs)
        @inbounds A[subs[i]] += val[i]
    end
    A
end

## Methods and Types Definitions


#1) Finds AKM largest connected set
"""
$(SIGNATURES)

Returns a tuple of observation belonging to the largest connected set with the corresponding identifiers and outcomes. This requires to have the data sorted by first identifier, and time period (e.g. we sort by worked id and year). This is also the set where we can run AKM models with the original data.

### Arguments
* `y`: outcome (e.g. log wage)
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
* `settings`: settings based on data type `VCHDFESettings`. Please see the reference provided below.
"""
function find_connected_set(y, first_idvar, second_idvar, settings)

    #Observation identifier to join later the FE
    obs_id = collect(1:size(y,1))

    seconds = unique(sort(second_idvar))
    second_id = indexin(second_idvar, seconds)

    firsts = unique(sort(first_idvar))
    first_id = indexin(first_idvar, firsts)

    #Save ids just in case
    second_id_old=second_idvar;
    first_id_old=first_idvar;

    n_seconds=length(seconds)
    n_firsts=length(firsts)
    second_id=second_id.+n_firsts

    graph_size=length(seconds)+length(firsts)
    G=Graph(graph_size)
    for i in 1:size(y,1)
        add_edge!(G, second_id[i], first_id[i])
    end

    cs=connected_components(G)
    (settings.print_level > 1) && println("Largest connected set has been found")
    pos = indexin(  [  maximum(size.(cs,1))] , size.(cs,1) )[1]
    lcs=cs[pos]

    connected_seconds=lcs[lcs.>n_firsts]

    sel=findall(in(connected_seconds),second_id)

    obs_id = [obs_id[x] for x in sel ]
    yvec = [y[x] for x in sel ]
    second_id = [second_idvar[x] for x in sel ]
    first_id = [first_idvar[x] for x in sel ]

    #Relabel seconds
    seconds = unique(sort(second_id))
    second_id = indexin(second_id, seconds)

    #Relabel firsts
    first_ids = unique(sort(first_id))
    first_id = indexin(first_id, first_ids)

    return (obs_id = obs_id , y = yvec , first_id = first_id, second_id = second_id )

end

#2) Pruning and finding Leave-Out Largest connected set
"""
$(SIGNATURES)

This function prunes the dataset from articulation points. If the first identifier is worker id it means that it prunes workers that would disconnect the graph if they were dropped.

### Arguments
* `yvec`: outcome (e.g. log wage)
* `first_idvar`: first identifier (e.g. worker id)
* `second_idvar`: second identifier (e.g. firm id)
* `obs_id`: observation identifier.
* `settings`: settings based on data type `VCHDFESettings`. Please see the reference provided below.
"""
function prunning_connected_set(yvec, first_idvar, second_idvar, obs_id, settings)

    seconds = unique(sort(second_idvar))
    second_id = indexin(second_idvar, seconds)
    firsts = unique(sort(first_idvar))
    first_id = indexin(first_idvar, firsts)

    nbadfirsts=1
    while nbadfirsts>0

        n_seconds=length(seconds)
        n_firsts=length(firsts)
        second_id=second_id.+n_firsts

        graph_size=n_firsts + n_seconds
        lcs_graph=Graph(graph_size)
            for i in 1:size(yvec,1)
                add_edge!(lcs_graph, second_id[i], first_id[i])
            end

        #get articulation vertex
        artic_vertex = articulation(lcs_graph)

        sel=findall(x->x==nothing, indexin(first_id,artic_vertex))

        bad_firsts = artic_vertex[artic_vertex.<=n_firsts]
        nbadfirsts = size(bad_firsts,1)

        settings.print_level > 1 && println("Number of $(settings.first_id_display_small) that are articulation points: ", nbadfirsts)

        #Restrict the sample
        yvec = [yvec[x] for x in sel ]
        second_id = [second_id[x] for x in sel ]
        first_id = [first_id[x] for x in sel ]
        obs_id = [obs_id[x] for x in sel ]

        #Relabel seconds
        seconds = unique(sort(second_id))
        second_id = indexin(second_id, seconds)

        #Relabel firsts
        first_ids = unique(sort(first_id))
        first_id = indexin(first_id, first_ids)

        n_firsts=maximum(first_id)
        n_seconds=maximum(second_id);

        second_id = second_id .+ n_firsts

        #Constructing new Graph
        G=Graph(n_firsts+n_seconds)
        Nprimeprime = size(sel,1)
        for i in 1:Nprimeprime
            add_edge!(G, second_id[i], first_id[i])
        end

        #Find Largest Connected Set (LCS)
        cs=connected_components(G)

        pos = indexin(  [  maximum(size.(cs,1))] , size.(cs,1) )[1]
        lcs=cs[pos]

        connected_seconds=lcs[lcs.>n_firsts]


        sel=findall(in(connected_seconds),second_id)

        obs_id = [obs_id[x] for x in sel ]
        yvec = [yvec[x] for x in sel ]
        second_id = [second_id[x] for x in sel ]
        first_id = [first_id[x] for x in sel ]

        #Relabel seconds
        seconds = unique(sort(second_id))
        second_id = indexin(second_id, seconds)

        #Relabel firsts
        first_ids = unique(sort(first_id))
        first_id = indexin(first_id, first_ids)

    end

    return (obs_id = obs_id , y = yvec , first_id = first_id, second_id = second_id )
end

#3) Drops Single Observations
"""
$(SIGNATURES)

This function drops observations that correspond with first identifiers with a single observation. For example, if first identifier is worker id, it will drop observations for workers that only appear once in the data.

### Arguments
* `yvec`: outcome (e.g. log wage)
* `first_idvar`: first identifier (e.g. worker id)
* `second_idvar`: second identifier (e.g. firm id)
* `obs_id`: observation identifier.
"""
function drop_single_obs(yvec, first_idvar, second_idvar,obs_id)

    seconds = unique(sort(second_idvar))
    second_id = indexin(second_idvar, seconds)
    firsts = unique(sort(first_idvar))
    first_id = indexin(first_idvar, firsts)

    T=accumarray(first_id,ones(length(first_id)))
    T = [T[x] for x in first_id]
    sel = T.>1
    sel = findall(x->x==true,sel )

    obs_id = [obs_id[x] for x in sel ]
    yvec = [yvec[x] for x in sel ]
    second_id = [second_id[x] for x in sel ]
    first_id = [first_id[x] for x in sel ]

    #Relabel seconds
    seconds = unique(sort(second_id))
    second_id = indexin(second_id, seconds)

    #Relabel firsts
    first_ids = unique(sort(first_id))
    first_id = indexin(first_id, first_ids)

    return (obs_id = obs_id , y = yvec , first_id = first_id, second_id = second_id )
end

#4) Finds the observations connected for every cluster
function index_constr(clustering_var, id, match_id )
    NT = length(clustering_var)
    counter = ones(size(clustering_var,1));

    #Indexing obs number per worked/id
    gcs = Int.(@transform(groupby(DataFrame(counter = counter, id = id), :id), gcs = cumsum(:counter)).gcs);
    maxD = maximum(gcs);

    index = collect(1:NT);

    #This will be the output, and we will append observations to it in the loop
    list_final=DataFrame(row = Int64[],col = Int64[], match_id = Int64[], id_cluster = Int64[]);


    for t=1:maxD
    rowsel =  findall(x->x==true,gcs.==t);
    list_base = DataFrame( id_cluster= [clustering_var[x] for x in rowsel], row =
    [index[x] for x in rowsel] , match_id = [match_id[x] for x in rowsel]  );

        for tt=t:maxD
            colsel =  findall(x->x==true,gcs.==tt);
            list_sel =  DataFrame( id_cluster= [clustering_var[x] for x in colsel], col =
            [index[x] for x in colsel] );

            merge = outerjoin(list_base, list_sel, on = :id_cluster)
            merge = dropmissing(merge)
            merge = merge[:,[:row, :col, :match_id, :id_cluster]]

            append!(list_final, merge)

        end
    end

    sort!(list_final, (:row));

    return Matrix(list_final)

end

#5) Compute Movers
"""
$(SIGNATURES)

Returns a vector that indicates whether the `first_id` (e.g. worker) is a mover across `second_id` (e.g. firms), as well as a vector with the number of periods that each `first_id` appears.

### Arguments
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
"""
function compute_movers(first_id,second_id)

    gcs = [NaN; first_id[1:end-1]]
    gcs = first_id.!=gcs

    lagsecond_id=[NaN; second_id[1:end-1]]
    for x in findall(x->x ==true , gcs)
        lagsecond_id[x] = NaN
    end

    stayer=(second_id.==lagsecond_id)
    for x in findall(x->x ==true , gcs)
        stayer[x] = true
    end

    stayer=Int.(stayer)
    stayer=accumarray(first_id,stayer)
    T=accumarray(first_id,ones(length(first_id)))
    stayer=T.==stayer
    movers=stayer.==false

    movers = [movers[x] for x in first_id]
    T = [T[x] for x in first_id]


    return (movers = movers, T = T)
end


#6) Eff res : Compute Effective Resistance - Lambda Matrices
"""
$(SIGNATURES)

This function computes the diagonal matrices containing Pii and Bii under the Exact Algorithm. See appendix in KSS for more information.

### Arguments
* `X`: the design matrix in the linear model.
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
* `match_id`: match identifier between first and second identifier.For example, match can be the identifier of every worker-firm combination.
* `K`: number of covariates in addition to the fixed effects. Currently only 0 is supported.
* `settings`: settings based on data type `VCHDFESettings`. Please see the reference provided below.
"""
function eff_res(::ExactAlgorithm, X,first_id,second_id,match_id, K, settings)

    #Indexing Observations
    elist = index_constr(collect(1:length(first_id)), first_id, match_id )

    #Dimensions
    NT = size(X,1)
    M = size(elist,1)
    J = maximum(second_id)
    N = maximum(first_id)

    #Define solver
    S_xx = X'*X

    # Create the solvers
    ldli, la = computeLDLinv(S_xx)
    buffs = zeros(size(la)[1],Threads.nthreads())
    compute_sol = []
    for i in 1:Threads.nthreads()
        P = approxcholOperator(ldli,buffs[:,i])
        push!(compute_sol,approxcholSolver(P,la))
    end

    #Initialize output
    Pii = zeros(M)
    Bii_second = zeros(M)
    Bii_cov= settings.cov_effects ==true ? zeros(M) : nothing
    Bii_first= settings.first_id_effects == true ? zeros(M) : nothing

    #No controls case: We compute Pii,Bii for stayers manually
    if K == 0

        #Compute Auxiliaries
        movers , T = compute_movers(first_id, second_id)

        Nmatches = maximum(match_id)
        match_id_movers = [match_id[x] for x in findall(x->x==true, movers)]
        second_id_movers = [second_id[x] for x in findall(x->x==true, movers)]
        first_id_movers = [first_id[x] for x in findall(x->x==true, movers)]

        sel = unique(z -> match_id_movers[z], 1:length(match_id_movers))
        match_id_movers = match_id_movers[sel]
        second_id_movers = second_id_movers[sel]
        first_id_movers = first_id_movers[sel]

        maxT = maximum([T[x] for x in findall(x->x == false, movers)])

        counter = ones(Int,NT)
        #not sure if I can change id to first_id (in the arguments)
        gcs = Int.(@transform(groupby(DataFrame(counter = counter, first_id = first_id), :first_id), gcs = cumsum(:counter)).gcs)
        sel_stayers = (gcs.==1).*(movers.==false)
        stayers_matches_sel = [match_id[z] for z in findall(x->x == true , sel_stayers)]
        Tinv = 1 ./T
        elist_JLL = [first_id_movers N.+second_id_movers first_id_movers N.+second_id_movers]

        M = size(elist_JLL,1)
        Pii_movers = zeros(M)
        Bii_second_movers = zeros(M)
        Bii_cov_movers= settings.cov_effects ==true ? zeros(M) : nothing
        Bii_first_movers= settings.first_id_effects == true ? zeros(M) : nothing

        #Initializing dependent variables for solver
        Xright = sparse(collect(1:M),elist_JLL[:,1],1.0,M,N+J)
        Xright = Xright .+ sparse(collect(1:M),elist_JLL[:,2],-1.0,M,N+J)
        # N+J x N+J-1 restriction matrix
        S= sparse(1.0I, J-1, J-1)
        S=vcat(S,sparse(-zeros(1,J-1)))

        Xright = hcat(Xright[:,1:N], Xright[:,N+1:end]*S)

        Threads.@threads for i=1:M

            #Only one inversion needed for exact alg
            zexact = compute_sol[Threads.threadid()]( [Xright[i,:]...] ; verbose=false)

            #Compute Pii
            Pii_movers[i] = Xright[i,:]'*zexact

            #Compute Bii for seconds
            aux_right = zexact[N+1:N+J-1]
            aux_left = zexact[N+1:N+J-1]

            COV = cov(X[:,N+1:N+J-1]*aux_left,X[:,N+1:N+J-1]*aux_right)
            Bii_second_movers[i] = COV[1]*(NT-1)

            if Bii_first != nothing
                aux_right = zexact[1:N]
                aux_left = zexact[1:N]
                COV = cov(X[:,1:N]*aux_left,X[:,1:N]*aux_right)
                Bii_first_movers[i] = COV[1]*(NT-1)
            end

            if Bii_cov != nothing
                aux_right = zexact[N+1:N+J-1]
                aux_left = zexact[1:N]
                COV = cov(X[:,1:N]*aux_left,X[:,N+1:N+J-1]*aux_right)
                Bii_cov_movers[i] = COV[1]*(NT-1)
            end

        end

        (settings.print_level > 1) && println("Pii and Bii have been computed for movers.")
        #Assign Step
        Pii_movers = sparse(match_id_movers,ones(Int,length(match_id_movers)),Pii_movers[:,1],Nmatches,1)
        Pii_stayers = sparse(stayers_matches_sel,ones(Int,length(stayers_matches_sel)),[Tinv[x] for x in findall(x->x==true,sel_stayers)],Nmatches,1)
        Pii = Pii_movers.+Pii_stayers

        Bii_second = sparse(match_id_movers,ones(Int,length(match_id_movers)),Bii_second_movers[:,1],Nmatches,1)

        if settings.cov_effects == true
            Bii_cov = sparse(match_id_movers,ones(Int,length(match_id_movers)),Bii_cov_movers[:,1],Nmatches,1)
        end

        if settings.first_id_effects == true
            Bii_first = sparse(match_id_movers,ones(Int,length(match_id_movers)),Bii_first_movers[:,1],Nmatches,1)
            stayers = .!movers

            Threads.@threads for t=2:maxT #T=1 have Pii=1 so need to be dropped.

                sel = (gcs.==true).*stayers.*(T.==t)
                N_sel = sum(sel)

                if N_sel > 0
                    index_sel = findall(x->x==true,sel)
                    match_sel_aux = Int.([match_id[z] for z in index_sel])
                    first = index_sel[1]
                    Xuse = X[first,:]

                    ztilde = compute_sol[Threads.threadid()]([Xuse...] ;verbose=false)

                    aux_right = ztilde[1:N]
                    aux_left = ztilde[1:N]

                    COV = cov(X[:,1:N]*aux_left,X[:,1:N]*aux_right)
                    Bii_first_stayers = COV[1]*(NT-1)

                    Bii_first_stayers = sparse(match_sel_aux,ones(Int,length(match_sel_aux)),Bii_first_stayers,Nmatches,1)
                    Bii_first = Bii_first.+Bii_first_stayers
                end
            end
        end


    #Controls case
    elseif K> 0
        #AUX = initialize_auxiliary_variables(settings.lls_algorithm, X, elist, M,NT, N, J, K, settings)
        nparameters = N + J + K

        D=sparse(collect(1:NT),first_id,1)
        F=sparse(collect(1:NT),second_id,1)
        S= sparse(1.0I, J-1, J-1)
        S=vcat(S,sparse(-zeros(1,J-1)))

        Dvar = hcat(  D, spzeros(NT,nparameters-N) )
        Fvar = hcat(spzeros(NT,N), -F*S, spzeros(NT,nparameters-N-J) )
        #Wvar = hcat(spzeros(NT,N+J), controls )
        Xleft = X[elist[:,1],:]
        Xright = X[elist[:,2],:]

        Threads.@threads for i=1:M

                #Again, one inversion needed
                zexact = compute_sol[Threads.threadid()]([Xright[i,:]...];verbose=false)

                #Compute Pii
                Pii[i] = Xleft[i,:]'*zexact

                #Compute Bii for seconds
                aux_right = zexact[N+1:N+J-1,:]
                aux_left = zexact[N+1:N+J-1,:]

                COV = cov(X[:,N+1:N+J-1]*aux_left,X[:,N+1:N+J-1]*aux_right)
                Bii_second[i] = COV[1]*(NT-1)

                if Bii_first != nothing
                    aux_right = zexact[1:N]
                    aux_left = zexact[1:N]
                    COV = cov(X[:,1:N]*aux_left,X[:,1:N]*aux_right)
                    Bii_first[i] = COV[1]*(NT-1)
                end

                if Bii_cov != nothing
                    aux_right = zexact[N+1:N+J-1]
                    aux_left = zexact[1:N]
                    COV = cov(X[:,1:N]*aux_left,X[:,N+1:N+J-1]*aux_right)
                    Bii_cov[i] = COV[1]*(NT-1)
                end

        end

    end

    #Create matrices
    rows = elist[:,1]
    cols = elist[:,2]
    index_cluster = match_id

    #Censor
    Pii[ findall(Pii.>=0.99)] .= 0.99

    if K==0
        Pii = [Pii[x] for x in index_cluster]
        Bii_second = [Bii_second[x] for x in index_cluster]

        if settings.cov_effects == true
            Bii_cov = [Bii_cov[x] for x in index_cluster]
        end

        if settings.first_id_effects == true
            Bii_first = [Bii_first[x] for x in index_cluster]
        end

    end


    #Lambda P
    Lambda_P=sparse(rows,cols,Pii,NT,NT)
    Lambda_P=Lambda_P+triu(Lambda_P,1)'

    #Lambda B var(fe)
    Lambda_B_second=sparse(rows,cols,Bii_second,NT,NT)
    Lambda_B_second=Lambda_B_second+triu(Lambda_B_second,1)'

    Lambda_B_cov = nothing
    #Lambda B cov(fe,pe)
    if settings.cov_effects == true
        Lambda_B_cov=sparse(rows,cols,Bii_cov,NT,NT)
        Lambda_B_cov=Lambda_B_cov+triu(Lambda_B_cov,1)'
    end
    
    Lambda_B_first = nothing
    #Lambda B, var(pe)
    if settings.first_id_effects == true
        Lambda_B_first=sparse(rows,cols,Bii_first,NT,NT)
        Lambda_B_first=Lambda_B_first+triu(Lambda_B_first,1)'
    end


    #TODO: maybe we can make the function to be inplace with those Lambdas
    # if settings.first_id_effects == false & settings.cov_effects == false
    #     return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second)
    # elseif settings.first_id_effects == true & settings.cov_effects == false
    #     return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first)
    # elseif settings.first_id_effects == true  & settings.cov_effects == true
    #     return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first, Lambda_B_cov=Lambda_B_cov)
    # end
    return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first, Lambda_B_cov=Lambda_B_cov)
end


"""
$(SIGNATURES)

This function computes the diagonal matrices containing Pii and Bii under Johnson-Linderstrauss Algorithm. See appendix in KSS for more information.

### Arguments
* `X`: the design matrix in the linear model.
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
* `match_id`: match identifier between first and second identifier.For example, match can be the identifier of every worker-firm combination.
* `K`: number of covariates in addition to the fixed effects. Currently only 0 is supported.
* `settings`: settings based on data type `VCHDFESettings`. Please see the reference provided below.
"""
function eff_res(lev::JLAAlgorithm, X,first_id,second_id,match_id, K, settings)

    #Indexing Observations
    elist = index_constr(collect(1:length(first_id)), first_id, match_id )

    #Dimensions
    NT=size(X,1)
    M=size(elist,1)
    J = maximum(second_id)
    N = maximum(first_id)
    p = lev.num_simulations == 0 ? ceil(log(N+J)/0.01) : lev.num_simulations

    #Define solver
    S_xx = X'*X
    ldli, la = computeLDLinv(S_xx)
    buffs = zeros(size(la)[1],Threads.nthreads())
    compute_sol = []
    for i in 1:Threads.nthreads()
        P = approxcholOperator(ldli,buffs[:,i])
        push!(compute_sol,approxcholSolver(P,la;tol=1e-12))
    end

    #Pre-allocate solution vectors
    Z = zeros(N+J,Threads.nthreads())
    ZB = zeros(N+J,Threads.nthreads())
    ZB_first = settings.cov_effects ==true ? zeros(N+J,Threads.nthreads()) : nothing 

    #Pre-allocate Rademacher Vectors
    rademach = zeros(Threads.nthreads(), NT)
    
    #Initialize output
    Pii=zeros(M)
    Bii_second=zeros(M)
    Bii_cov= settings.cov_effects ==true ? zeros(M) : nothing
    Bii_first= settings.first_id_effects == true ? zeros(M) : nothing

    #No controls case
    if K == 0

        #Compute Auxiliaries
        movers , T = compute_movers(first_id, second_id)

        Nmatches = maximum(match_id)
        match_id_movers = [match_id[x] for x in findall(x->x==true, movers)]
        second_id_movers = [second_id[x] for x in findall(x->x==true, movers)]
        first_id_movers = [first_id[x] for x in findall(x->x==true, movers)]

        sel = unique(z -> match_id_movers[z], 1:length(match_id_movers))
        match_id_movers=match_id_movers[sel]
        second_id_movers=second_id_movers[sel]
        first_id_movers=first_id_movers[sel]

        maxT = maximum([T[x] for x in findall(x->x == false, movers)])

        counter = ones(Int,NT)
        #again, not sure about renaming the vars
        gcs = Int.(@transform(groupby(DataFrame(counter = counter, first_id = first_id), :first_id), gcs = cumsum(:counter)).gcs)
        sel_stayers=(gcs.==1).*(movers.==false)
        stayers_matches_sel=[match_id[z] for z in findall(x->x == true , sel_stayers)]
        Tinv=1 ./T
        elist_JLL=[first_id_movers N.+second_id_movers]

        M=size(elist_JLL,1)
        Pii_movers=zeros(M)
        Bii_second_movers=zeros(M)
        Bii_cov_movers= settings.cov_effects ==true ? zeros(M) : nothing
        Bii_first_movers= settings.first_id_effects == true ? zeros(M) : nothing

        #Initializing dependent variables for LSS
        Fvar= hcat(spzeros(NT,N), X[:,N+1:N+J-1])
        Dvar=hcat(X[:,1:N], spzeros(NT,J-1))

        (settings.print_level > 0) && println("Running JLA with ",p," simulations.")

        Threads.@threads for i=1:p
            #Draw Rademacher entry
            rademach[Threads.threadid(),:] =  rand(1,NT) .> 0.5
            rademach[Threads.threadid(),:]  = rademach[Threads.threadid(),:]  .- (rademach[Threads.threadid(),:]  .== 0)
            rademach[Threads.threadid(),:]  = rademach[Threads.threadid(),:]  ./sqrt(p)

            Z[1:end-1,Threads.threadid()]  .= compute_sol[Threads.threadid()]( [rademach[Threads.threadid(),:]'*X...] ; verbose=false)

            rademach[Threads.threadid(),:] = rademach[Threads.threadid(),:] .- mean(rademach[Threads.threadid(),:])
            ZB[1:end-1,Threads.threadid()] .= compute_sol[Threads.threadid()]( [rademach[Threads.threadid(),:]'*Fvar...] ; verbose=false)

            if settings.first_id_effects == true | settings.cov_effects == true
                ZB_first[1:end-1,Threads.threadid()] .= compute_sol[Threads.threadid()]( [rademach[Threads.threadid(),:]'*Dvar...] ; verbose=false)
            end

            #Z[:,Threads.threadid()] = [Z[:,Threads.threadid()];0.0]
            #ZB[:,Threads.threadid()] = [ZB[:,Threads.threadid()];0.0]
            #ZB_first[:,Threads.threadid()] = [ZB_first[:,Threads.threadid()];0.0]

            #Computing
            Pii_movers = Pii_movers .+ ( [Z[j,Threads.threadid()]  for j in elist_JLL[:,1] ]  .- [Z[j,Threads.threadid()]  for j in elist_JLL[:,2] ] ) .* ( [Z[j,Threads.threadid()]  for j in elist_JLL[:,1] ]  .- [Z[j,Threads.threadid()]  for j in elist_JLL[:,2] ] )
            Bii_second_movers = Bii_second_movers .+ ( [ZB[j,Threads.threadid()]  for j in elist_JLL[:,1] ]  .- [ZB[j,Threads.threadid()]  for j in elist_JLL[:,2] ] ) .* ( [ZB[j,Threads.threadid()]  for j in elist_JLL[:,1] ]  .- [ZB[j,Threads.threadid()]  for j in elist_JLL[:,2] ] )

            if settings.first_id_effects == true
                Bii_first_movers = Bii_first_movers .+  ( [ZB_first[j,Threads.threadid()] for j in elist_JLL[:,1] ]  .- [ZB_first[j,Threads.threadid()]  for j in elist_JLL[:,2] ] ) .* ( [ZB_first[j,Threads.threadid()]  for j in elist_JLL[:,1] ]  .- [ZB_first[j,Threads.threadid()] for j in elist_JLL[:,2] ] )
            end

            if settings.cov_effects == true
                Bii_cov_movers = Bii_cov_movers .+ ( [ZB[j,Threads.threadid()]   for j in elist_JLL[:,1] ]  .- [ZB[j,Threads.threadid()]  for j in elist_JLL[:,2] ] ) .* ( [ZB_first[j,Threads.threadid()]  for j in elist_JLL[:,1] ]  .- [ZB_first[j,Threads.threadid()]  for j in elist_JLL[:,2] ] )
            end

        end

        (settings.print_level > 1) && println("Pii and Bii have been computed for movers.")
        #Assign Step
        Pii_movers=sparse(match_id_movers,ones(Int,length(match_id_movers)),Pii_movers[:,1],Nmatches,1)
        Pii_stayers=sparse(stayers_matches_sel,ones(Int,length(stayers_matches_sel)),[Tinv[x] for x in findall(x->x==true,sel_stayers)],Nmatches,1)
        Pii=Pii_movers.+Pii_stayers

        Bii_second=sparse(match_id_movers,ones(Int,length(match_id_movers)),Bii_second_movers[:,1],Nmatches,1)

        if settings.cov_effects == true
            Bii_cov=sparse(match_id_movers,ones(Int,length(match_id_movers)),Bii_cov_movers[:,1],Nmatches,1)
        end

        if settings.first_id_effects == true
            Bii_first=sparse(match_id_movers,ones(Int,length(match_id_movers)),Bii_first_movers[:,1],Nmatches,1)
            stayers = .!movers

            Threads.@threads for t=2:maxT #T=1 have Pii=1 so need to be dropped.
                sel=(gcs.==true).*stayers.*(T.==t)
                N_sel=sum(sel)
                if N_sel > 0
                    index_sel=findall(x->x==true,sel)
                    match_sel_aux=Int.([match_id[z] for z in index_sel])
                    first=index_sel[1]
                    Xuse=X[first,:]

                    ztilde = compute_sol[Threads.threadid()]([X[first,:]...] ; verbose=false)

                    aux_right=ztilde[1:N]
                    aux_left=ztilde[1:N]

                    COV=cov(X[:,1:N]*aux_left,X[:,1:N]*aux_right)
                    Bii_first_stayers=COV[1]*(NT-1)

                    Bii_first_stayers=sparse(match_sel_aux,ones(Int,length(match_sel_aux)),Bii_first_stayers,Nmatches,1)
                    Bii_first=Bii_first.+Bii_first_stayers
                end
            end
        end



    elseif K> 0
        #AUX = initialize_auxiliary_variables(settings.lls_algorithm, X, elist, M,NT, N, J, K, settings)
        nparameters = N + J + K

        D=sparse(collect(1:NT),first_id,1)
        F=sparse(collect(1:NT),second_id,1)
        Dvar = hcat(  D, spzeros(NT,nparameters-N) )
        Fvar = hcat(spzeros(NT,N), F, spzeros(NT,nparameters-N-J) )
        #Wvar = hcat(spzeros(NT,N+J), controls )
        Xleft = X[elist[:,1],:]
        Xright = X[elist[:,2],:]

        Threads.@threads for i=1:p

            #Rademacher Entries
            rademach = rand(1,NT) .> 0.5
            rademach = rademach - (rademach .== 0)
            rademach = rademach ./sqrt(p)


            Zleft = compute_sol[Threads.threadid()]( [rademach*Xleft...] ; verbose=false)
            #Zright = lss(settings.lls_algorithm, X, rademach*Xright, settings)

            Pii = Pii .+ (X*Zleft).^2

            rademach = rademach .- mean(rademach)

            aux = compute_sol[Threads.threadid()]( [rademach*Fvar...] ;verbose=false)
            ZF = X*aux

            Bii_second = Bii_second .+ ZF.^2 ./NT

            if settings.first_id_effects == true |    settings.cov_effects == true
                aux = compute_sol[Threads.threadid()]( [rademach*Dvar...] ;verbose=false)
                ZD = X*aux
            end

            if settings.first_id_effects==true
                Bii_first = Bii_first .+ (ZD).^2 ./ NT
            end

            if settings.cov_effects==true
                Bii_cov = Bii_cov .+ (ZD .* ZF ) ./ NT
            end

        end


    end

    #Create matrices
    rows = elist[:,1]
    cols = elist[:,2]
    index_cluster = elist[:,3]

    #Censor
    Pii[ findall(Pii.>=0.99)] .= 0.99

    if K==0
        Pii = [Pii[x] for x in index_cluster]
        Bii_second = [Bii_second[x] for x in index_cluster]

        if settings.cov_effects == true
            Bii_cov = [Bii_cov[x] for x in index_cluster]
        end

        if settings.first_id_effects == true
            Bii_first = [Bii_first[x] for x in index_cluster]
        end

    end


    #Lambda P
    Lambda_P=sparse(rows,cols,Pii,NT,NT)
    Lambda_P=Lambda_P+triu(Lambda_P,1)'

    #Lambda B var(second effects)
    Lambda_B_second=sparse(rows,cols,Bii_second,NT,NT)
    Lambda_B_second=Lambda_B_second+triu(Lambda_B_second,1)'

    Lambda_B_cov = nothing
    #Lambda B cov(fe,pe)
    if settings.cov_effects == true
        Lambda_B_cov=sparse(rows,cols,Bii_cov,NT,NT)
        Lambda_B_cov=Lambda_B_cov+triu(Lambda_B_cov,1)'
    end

    #Lambda B, var(pe)
    Lambda_B_first = nothing
    if settings.first_id_effects == true
        Lambda_B_first=sparse(rows,cols,Bii_first,NT,NT)
        Lambda_B_first=Lambda_B_first+triu(Lambda_B_first,1)'
    end


    # if settings.first_id_effects == false & settings.cov_effects == false
    #     return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second)
    # elseif settings.first_id_effects == true & settings.cov_effects == false
    #     return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first)
    # elseif settings.first_id_effects == true  & settings.cov_effects == true
    #     return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first, Lambda_B_cov=Lambda_B_cov)
    # end

    return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first, Lambda_B_cov=Lambda_B_cov)

end



#9) Creates match first_id using second_id id
"""
$(SIGNATURES)

Computes a match identifier for every combination of first and second identifier. For example, this can be the match identifier of worker-firm combinations.

### Arguments
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
"""
function compute_matchid(second_id,first_id)
    match_id2 = string.(first_id).*string.("+",second_id)
    match_id2 = indexin(match_id2, unique(match_id2))

    return match_id2
end

"""
$(SIGNATURES)

Returns the bias-corrected components, the vector of coefficients, the corresponding fixed effects for every observation, and the diagonal matrices containing the Pii and Biis. 

### Arguments
* `y`: outcome vector
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
* `settings`: settings based on `VCHDFESettings`
* `controls`: at this version only `controls=nothing` is supported.
"""
function leave_out_estimation(y,first_id,second_id,controls,settings)

    #Create matrices for computations
    NT = size(y,1)
    J = maximum(second_id)
    N = maximum(first_id)
    K = controls ==nothing ? 0 : size(controls,2)
    nparameters = N + J + K

    match = compute_matchid(first_id, second_id)

    #first (worker) Dummies
    D = sparse(collect(1:NT),first_id,1)

    #second (firm) Dummies
    F = sparse(collect(1:NT),second_id,1)

    # N+J x N+J-1 restriction matrix
    S= sparse(1.0I, J-1, J-1)
    S=vcat(S,sparse(-zeros(1,J-1)))

    X = hcat(D, -F*S)

    #SET DIMENSIONS
    n=size(X,1)

    # PART 1: ESTIMATE HIGH DIMENSIONAL MODEL
    xx=X'*X
    xy=X'*y
    compute_sol = approxcholSolver(xx;verbose = settings.print_level > 0)
    beta = compute_sol([xy...];verbose=false)
    eta=y-X*beta

    pe=D * beta[1:N]
    fe=F*S * beta[N+1:N+J-1]

    println("\n","Plug-in Variance Components:")

    σ2_ψ_AKM = var(fe)
    println("Plug-in Variance of $(settings.second_id_display_small) Effects: ", σ2_ψ_AKM )
    σ2_α_AKM = var(pe)
    println("Plug-in Variance of $(settings.first_id_display_small) Effects: ", σ2_α_AKM )
    σ2_ψα_AKM = cov(pe,-fe)
    println("Plug-in Covariance of $(settings.first_id_display_small)-$(settings.second_id_display_small) Effects: ", σ2_ψα_AKM, "\n")


    #Part 2: Compute Pii, Bii
    @unpack Lambda_P, Lambda_B_second,Lambda_B_first, Lambda_B_cov = eff_res(settings.leverage_algorithm, X,first_id,second_id,match, K, settings)

    #Compute Leaveo-out residual
    I_Lambda_P = I-Lambda_P
    eta_h = I_Lambda_P\eta

    #Compute bias corrected variance comp of second (Firm) Effects
    θ_second = σ2_ψ_AKM -(1/NT)*y'*Lambda_B_second*eta_h

    θ_first = settings.first_id_effects==true  ? σ2_α_AKM -(1/NT)*y'*Lambda_B_first*eta_h : nothing

    θCOV = settings.cov_effects==true ? σ2_ψα_AKM -(1/NT)*y'*Lambda_B_cov*eta_h : nothing

    Pii = diag(Lambda_P)

    Bii_second = diag(Lambda_B_second)

    Bii_first = settings.first_id_effects == true ? diag(Lambda_B_first) : nothing

    Bii_cov = settings.cov_effects == true ? diag(Lambda_B_cov) : nothing

    #TODO print estimates
    if settings.print_level > 0
        println("Bias-Corrected Variance Components:")
        println("Bias-Corrected variance of $(settings.second_id_display_small): ", θ_second)
        (settings.first_id_effects > 0) && println("Bias-Corrected variance of $(settings.first_id_display_small): ", θ_first)
        (settings.cov_effects > 0) && println("Bias-Corrected covariance of $(settings.first_id_display_small)-$(settings.second_id_display_small) effects: ", θCOV)
    end

    return (θ_first = θ_first, θ_second = θ_second, θCOV = θCOV, β = beta, Dalpha = pe, Fpsi = fe, Pii = Pii, Bii_first = Bii_first,
            Bii_second = Bii_second, Bii_cov = Bii_cov)

end


"""
$(SIGNATURES)

Returns a tuple with the observation number of the original dataset that belongs to the Leave-out connected set as described in Kline,Saggio, Solvesten. It also provides the corresponding outcome and identifiers in this connected set. 

### Arguments
* `y`: outcome vector
* `first_id`: first identifier (e.g. worker id)
* `second_id`: second identifier (e.g. firm id)
* `settings`: settings based on `VCHDFESettings`
* `controls`: at this version only `controls=nothing` is supported.
"""
function get_leave_one_out_set(y, first_id, second_id, settings, controls)
    # @assert settings.first_id_effects == true && settings.cov_effects == true

    # compute y, id firmid, controls, settings
    # compute y, first_id second_id, controls, settings
    (settings.print_level > 0) && println("Finding the leave-one-out connected set")
    @unpack obs_id,  y  , first_id , second_id  = find_connected_set(y,first_id,second_id,settings)
    @unpack obs_id,  y  , first_id , second_id  = prunning_connected_set(y,first_id,second_id, obs_id,settings)
    @unpack obs_id,  y  , first_id , second_id  = drop_single_obs(y,first_id,second_id, obs_id)
    controls == nothing ? nothing : controls = controls[obs_id,:]

    return (obs = obs_id, y = y, first_id = first_id, second_id = second_id, controls = controls)
end

function leave_out_KSS(y,first_id,second_id,controls,settings)

    @unpack obs,  y  , first_id , second_id, controls = get_leave_one_out_set(y, first_id, second_id, settings, nothing)

    @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation2(y,first_id,second_id,controls,settings)

end


function sigma_for_stayers(y,first_id, second_id, weight, b)
    #Go back to person-year space 
    first_id_weighted = vcat(fill.(first_id, weight)...)
    second_id_weighted = vcat(fill.(second_id, weight)...)

    #Compute Pii for Stayers = inverse of # of obs 
    T = accumarray(id, 1)
    Pii = 1 ./T
    Mii = 1 - Pii 

    #Compute OLS residual 
    NT = size(y,1)
    D = sparse(collect(1:NT),first_id_weighted, 1)
    F = sparse(collect(1:NT),second_id_weighted, 1)
    J = size(F,2)
    S= sparse(1.0I, J-1, J-1)
    S=vcat(S,sparse(-zeros(1,J-1)))
    X = hcat(D, -F*S)

    eta = y - X*b 
    eta_h = eta./Mii 
    sigma_stayers = (y - mean(y)).*eta_h

    #Collapse to match
    match_id = compute_matchid(second_id_weighted,first_id_weighted)
    sigma_stayers = (combine(groupby(DataFrame(sigma = sigma_stayers , match_id = match_id), :match_id), :sigma => mean).sigma_mean)
    
    return sigma_stayers
end


function kss_quadratic_form(sigma_i, A_1, A_2, beta, Bii)
    right                               = A_2*beta
    left                                = A_1*beta
    theta                               = cov(left,right)
    theta                               = theta[1]
    dof                                 = size(left,1)-1
    theta_KSS                           = theta-(1/dof)*sum(Bii.*sigma_i)
end

function leave_out_estimation2(y,first_id,second_id,controls,settings)

    #Create matrices for computations
    NT = size(y,1)
    J = maximum(second_id)
    N = maximum(first_id)
    K = controls ==nothing ? 0 : size(controls,2)
    nparameters = N + J + K

    match_id = compute_matchid(first_id, second_id)

    #first (worker) Dummies
    D = sparse(collect(1:NT),first_id,1)

    #second (firm) Dummies
    F = sparse(collect(1:NT),second_id,1)

    # N+J x N+J-1 restriction matrix
    S= sparse(1.0I, J-1, J-1)
    S=vcat(S,sparse(-zeros(1,J-1)))

    X = hcat(D, -F*S)

    #SET DIMENSIONS
    n=size(X,1)

    # PART 1: ESTIMATE HIGH DIMENSIONAL MODEL
    xx=X'*X
    xy=X'*y
    compute_sol = approxcholSolver(xx;verbose = settings.print_level > 0)
    beta = compute_sol([xy...];verbose=false)
    eta=y-X*beta

    pe=D * beta[1:N]
    fe=F*S * beta[N+1:N+J-1]

    println("\n","Plug-in Variance Components:")

    σ2_ψ_AKM = var(fe)
    println("Plug-in Variance of $(settings.second_id_display_small) Effects: ", σ2_ψ_AKM )
    σ2_α_AKM = var(pe)
    println("Plug-in Variance of $(settings.first_id_display_small) Effects: ", σ2_α_AKM )
    σ2_ψα_AKM = cov(pe,-fe)
    println("Plug-in Covariance of $(settings.first_id_display_small)-$(settings.second_id_display_small) Effects: ", σ2_ψα_AKM, "\n")

    #Part 2: Collapse & Reweight (if needed)
    weight = ones(NT,1) 
    y_py = y
    if settings.leave_out_level == "match"
        #Compute Weights 
        weight = accumarray(match_id, 1)

        first_id_weighted = Int.(transform(groupby(DataFrame(first_id = id, match_id = match_id), :match_id), :first_id => mean  => :first_id_weighted).first_id_weighted)
        second_id_weighted = Int.(transform(groupby(DataFrame(first_id = id, match_id = match_id), :match_id), :first_id => mean  => :first_id_weighted).first_id_weighted)

        #Collapse at match level averages
        first_id = Int.(combine(groupby(DataFrame(first_id = first_id, match_id = match_id), :match_id), :first_id => mean).first_id_mean)
        second_id = Int.(combine(groupby(DataFrame(second_id = second_id, match_id = match_id), :match_id), :second_id => mean).second_id_mean)
        y = (combine(groupby(DataFrame(y = y, match_id = match_id), :match_id), :y => mean).y_mean)

        #At person-year level (like a Matlab's repelem)
        first_id_weighted = vcat(fill.(first_id,weight)...)
        second_id_weighted = vcat(fill.(second_id,weight)...)

        #Build Design
        NT = size(y,1)
        J = maximum(second_id)
        N = maximum(first_id)
        nparameters = N + J + K
        D = sparse(collect(1:NT),first_id,1)
        F = sparse(collect(1:NT),second_id,1)
        S= sparse(1.0I, J-1, J-1)
        S=vcat(S,sparse(-zeros(1,J-1)))

        X = hcat(D, -F*S)   
        Dvar = hcat( sparse(collect(1:length(match_id)),first_id_weighted,1) , spzeros(NT,J-1))
        Fvar = hcat(spzeros(NT,N), -1*sparse(collect(1:length(match_id)), second_id_weighted,1  )*S )
            
        #Weighting
        weight_mat = sparse(collect(1:NT), collect(1:NT), weight.^(0.5) , NT, NT )
        X = weight_mat * X 
        y = weight_mat * y 

    else 
        #Build Design
        NT = size(y,1)
        J = maximum(second_id)
        N = maximum(first_id)
        nparameters = N + J + K
        D = sparse(collect(1:NT),first_id,1)
        F = sparse(collect(1:NT),second_id,1)
        S= sparse(1.0I, J-1, J-1)
        S=vcat(S,sparse(-zeros(1,J-1)))

        X = hcat(D, -F*S)
        Dvar = hcat(  D, spzeros(NT,J-1) )
        Fvar = hcat(spzeros(NT,N), -F*S )

    end

    #Part 3: Compute Pii, Bii
    @unpack Pii , Mii  , correction_JLA , Bii_first , Bii_second , Bii_cov = leverages(settings.leverage_algorithm, X, Dvar, Fvar, settings)

    (settings.print_level > 1) && println("Pii and Bii have been computed.")

    #Compute Leave-out residual
    xx=X'*X
    xy=X'*y
    compute_sol = approxcholSolver(xx;verbose = settings.print_level > 0)
    beta = compute_sol([xy...];verbose=false)
    eta=y-X*beta

    eta_h = eta ./ Mii
    sigma_i = ( ( y .- mean(y) ) .* eta_h ) .* correction_JLA

    if settings.leave_out_level == "match"
        T = accumarray(id,1)
        stayers = T.==1 
        stayers = [stayers[j] for j in id]

        sigma_stayers = sigma_for_stayers(y_py, first_id, second_id, weight, beta)
        sigma[stayers] .= [sigma_stayers[j] for j in stayers]
    end


    #Compute bias corrected variance comp of second (Firm) Effects
    θ_second = kss_quadratic_form(sigma_i, Fvar, Fvar, beta, Bii_second)

    θ_first = settings.first_id_effects==true  ? kss_quadratic_form(sigma_i, Dvar, Dvar, beta, Bii_first) : nothing

    θCOV = settings.cov_effects==true ? kss_quadratic_form(sigma_i, Fvar, Dvar, beta, Bii_cov) : nothing

    Pii = diag(Pii)

    Bii_second = diag(Bii_second)

    Bii_first = settings.first_id_effects == true ? diag(Bii_first) : nothing

    Bii_cov = settings.cov_effects == true ? diag(Bii_cov) : nothing

    #TODO print estimates
    if settings.print_level > 0
        println("Bias-Corrected Variance Components:")
        println("Bias-Corrected variance of $(settings.second_id_display_small): ", θ_second)
        (settings.first_id_effects > 0) && println("Bias-Corrected variance of $(settings.first_id_display_small): ", θ_first)
        (settings.cov_effects > 0) && println("Bias-Corrected covariance of $(settings.first_id_display_small)-$(settings.second_id_display_small) effects: ", θCOV)
    end

    return (θ_first = θ_first, θ_second = θ_second, θCOV = θCOV, β = beta, Dalpha = pe, Fpsi = fe, Pii = Pii, Bii_first = Bii_first,
            Bii_second = Bii_second, Bii_cov = Bii_cov)
end



function leverages(::ExactAlgorithm, X,Dvar,Fvar, settings)

    M = size(X,1)

    #Define solver
    S_xx = X'*X

    # Create the solvers
    ldli, la = computeLDLinv(S_xx)
    buffs = zeros(size(la)[1],Threads.nthreads())
    compute_sol = []
    for i in 1:Threads.nthreads()
        P = approxcholOperator(ldli,buffs[:,i])
        push!(compute_sol,approxcholSolver(P,la))
    end

    #Initialize output
    Pii = zeros(M)
    Bii_second = zeros(M)
    Bii_cov= settings.cov_effects ==true ? zeros(M) : nothing
    Bii_first= settings.first_id_effects == true ? zeros(M) : nothing

    Threads.@threads for i=1:M

        #Only one inversion needed for exact alg
        zexact = compute_sol[Threads.threadid()]( [X[i,:]...] ; verbose=false)

        #Compute Pii
        Pii[i] = X[i,:]'*zexact

        #Compute Bii 
        COV = cov(Fvar*zexact,Fvar*zexact)
        Bii_second[i] = COV[1]*(size(Fvar,1)-1)

        if Bii_first != nothing
            COV = cov(Fvar*zexact,Dvar*zexact)
            Bii_first[i] = COV[1]*(size(Fvar,1)-1)
        end

        if Bii_cov != nothing
            COV = cov(Dvar*zexact,Dvar*zexact)
            Bii_cov[i] = COV[1]*(size(Dvar,1)-1)
        end

    end

    #Censor
    Pii[ findall(Pii.>=0.99)] .= 0.99

    correction_JLA = 1 
    Mii = 1 .- Pii

    return (Pii = Pii , Mii = Mii , correction_JLA = correction_JLA, Bii_first = Bii_first , Bii_second = Bii_second , Bii_cov = Bii_cov)
end



function leverages(lev::JLAAlgorithm, X,Dvar,Fvar, settings)

    NT = size(X,1)
    FE = size(X,2)

    p = lev.num_simulations == 0 ? ceil(log(FE)/0.01) : lev.num_simulations

    #Define solver
    S_xx = X'*X

    # Create the solvers
    ldli, la = computeLDLinv(S_xx)
    buffs = zeros(size(la)[1],Threads.nthreads())
    compute_sol = []
    for i in 1:Threads.nthreads()
        P = approxcholOperator(ldli,buffs[:,i])
        push!(compute_sol,approxcholSolver(P,la))
    end

    #Clear that memory
    ldli = nothing 
    la = nothing 
    buffs = nothing
    S_xx = nothing

    #Pre-allocate Rademacher Seeds
    mts = MersenneTwister.(1:Threads.nthreads())
    
    #Initialize output
    Pii=zeros(NT)
    Bii_second=zeros(NT)
    Bii_cov= settings.cov_effects ==true ? zeros(NT) : nothing
    Bii_first= settings.first_id_effects == true ? zeros(NT) : nothing

    Pii_sq = zeros(NT)
    Mii_sq = zeros(NT)
    Pii_Mii = zeros(NT)

    Threads.@threads for i=1:p

        rademach =  rand(mts[Threads.threadid()],1,NT) .> 0.5
        rademach  = rademach  .- (rademach .== 0)
        rademach  = rademach / sqrt(p)

        Z= compute_sol[Threads.threadid()]( [rademach*X...] ; verbose=false)
        Z = X*Z

        #Auxiliaries for Non-linear Correction

        aux = Z.^2 / p 
        Pii = Pii .+ aux 

        aux = Z.^4 / p 
        Pii_sq = Pii_sq .+ aux 

        aux			= ((rademach .- Z).^2)/p
		Mii			= Mii .+ aux    
		aux			= ((rademach .- Z).^4)/p
		Mii_sq		= Mii_sq .+ aux
		
        Pii_Mii		= Pii_Mii .+ ((Z.^2).*((rademach .- Z).^2) )/p 
        
        #Demeaned Rademacher
        rademach =  rand(mts[Threads.threadid()],1,size(Fvar,1)) .> 0.5
        rademach  = rademach  .- (rademach .== 0)
        rademach  = rademach /sqrt(p)
        rademach = rademach .- mean(rademach)

        Z .= compute_sol[Threads.threadid()]( [rademach*Fvar...] ; verbose=false)
        Z = X*Z
        Bii_second = Bii_second .+ (Z.*Z) 

        if settings.first_id_effects == true | settings.cov_effects == true
            Z_pe .= compute_sol[Threads.threadid()]( [rademach*Dvar...] ; verbose=false)
            Z_pe = X*Z_pe

            if settings.first_id_effects == true 
                Bii_first = Bii_first .+ (Z_pe.*Z_pe)
            end

            if settings.cov_effects == true 
                Bii_cov = Bii_cov .+ (Z_pe.*Z)
            end

        end    
    end

    #Censor
    Pii[ findall(Pii.>=0.99)] .= 0.99

    #Account for Non-linear Bias
    Pii = Pii ./ (Pii .+ Mii)
    Mii = 1 .- Pii 
    Vi = (1/p)*((Mii.^2).*Pii_sq+(Pii.^2).*Mii_sq-2*Mii.*Pii.*Pii_Mii)
    Bi = (1/p)*(Mii.*Pii_sq-Pii.*Mii_sq+2*(Mii-Pii).*Pii_Mii)
    correction_JLA = (1-Vi./(Mii.^2)+Bi./Mii)       

    return (Pii = Pii , Mii = Mii , correction_JLA = correction_JLA, Bii_first = Bii_first , Bii_second = Bii_second , Bii_cov = Bii_cov)
end


## Dependencies

#To install Laplacians you need to put in Julia's Prompt "add Laplacians#master"
#For some reason Pkg.add() doesn't work with this.
#Pkg.add("MATLAB")
#using Laplacians
#using Pkg
#include(string(Pkg.dir("Laplacians") , "/src/matlabSolvers.jl"))
getlagged(x) = [NaN; x[1:(end - 1)]]


#Defining types and structures
abstract type AbstractLLSAlgorithm end
abstract type AbstractGraphLLSAlgorithm <: AbstractLLSAlgorithm end  # i.e. requires graph laplacian
struct CMGPreconditionedLLS <: AbstractGraphLLSAlgorithm  end
struct AMGPreconditionedLLS <: AbstractGraphLLSAlgorithm  end
struct DirectLLS <: AbstractLLSAlgorithm  end  # no graph required

abstract type AbstractLeverageAlgorithm end
struct ExactAlgorithm <: AbstractLeverageAlgorithm  end
@with_kw struct JLAAlgorithm <: AbstractLeverageAlgorithm
    num_simulations::Int64 = 0
end

@with_kw struct VCHDFESettings{LeverageAlgorithm}
    cg_maxiter::Int64 = 300
    leverage_algorithm::LeverageAlgorithm = ExactAlgorithm()
    #clustering_level::String = "obs"
    first_id_effects::Bool = false
    cov_effects::Bool = false
    print_level::Int64 = 1
    first_id_display_small::String = "person"
    first_id_display::String = "Person"
    second_id_display_small::String = "firm"
    second_id_display::String = "Firm"
    observation_id_display_small::String = "wage"
    observation_id_display::String = "Wage"
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
    compute_movers(first_id,second_id)

Returns a vector that indicates whether the `first_id` (e.g. worker) is a mover across `second_id` (e.g. firms), as well as a vector with the number of periods that each `first_id' appears.
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

    #Lambda B cov(fe,pe)
    if settings.cov_effects == true
        Lambda_B_cov=sparse(rows,cols,Bii_cov,NT,NT)
        Lambda_B_cov=Lambda_B_cov+triu(Lambda_B_cov,1)'
    end

    #Lambda B, var(pe)
    if settings.first_id_effects == true
        Lambda_B_first=sparse(rows,cols,Bii_first,NT,NT)
        Lambda_B_first=Lambda_B_first+triu(Lambda_B_first,1)'
    end


    #TODO: maybe we can make the function to be inplace with those Lambdas
    if settings.first_id_effects == false & settings.cov_effects == false
        return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second)
    elseif settings.first_id_effects == true & settings.cov_effects == false
        return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first)
    elseif settings.first_id_effects == true  & settings.cov_effects == true
        return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first, Lambda_B_cov=Lambda_B_cov)
    end

end



function eff_res(lev::JLAAlgorithm, X,first_id,second_id,match_id, K, settings)

    #Indexing Observations
    elist = index_constr(collect(1:length(first_id)), first_id, match_id )

    #Dimensions
    NT=size(X,1)
    M=size(elist,1)
    J = maximum(second_id)
    N = maximum(first_id)
    p = lev.num_simulations == 0 ? ceil(log2(NT)/0.005) : lev.num_simulations

    #Define solver
    S_xx = X'*X
    ldli, la = computeLDLinv(S_xx)
    buffs = zeros(size(la)[1],Threads.nthreads())
    compute_sol = []
    for i in 1:Threads.nthreads()
        P = approxcholOperator(ldli,buffs[:,i])
        push!(compute_sol,approxcholSolver(P,la))
    end

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
        elist_JLL=[first_id_movers N.+second_id_movers first_id_movers N.+second_id_movers]

        M=size(elist_JLL,1)
        Pii_movers=zeros(M)
        Bii_second_movers=zeros(M)
        Bii_cov_movers= settings.cov_effects ==true ? zeros(M) : nothing
        Bii_first_movers= settings.first_id_effects == true ? zeros(M) : nothing

        #Initializing dependent variables for LSS
        Fvar= hcat(spzeros(NT,N), X[:,N+1:N+J-1])
        Dvar=hcat(X[:,1:N], spzeros(NT,J-1))

        Threads.@threads for i=1:p

            #Draw Rademacher entry
            rademach = rand(1,NT) .> 0.5
            rademach = rademach - (rademach .== 0)
            rademach = rademach ./sqrt(p)

            Z  = compute_sol[Threads.threadid()]( [rademach*X...] ; verbose=false)

            rademach = rademach .- mean(rademach)
            ZB = compute_sol[Threads.threadid()]( [rademach*Fvar...] ; verbose=false)


            if settings.first_id_effects == true | settings.cov_effects == true
                ZB_first = compute_sol[Threads.threadid()]( [rademach*Dvar...] ; verbose=false)
            end

            Z = [Z;0.0]
            ZB = [ZB;0.0]
            ZB_first = [ZB_first;0.0]

            #Computing
            Pii_movers = Pii_movers .+ ( [Z[j]  for j in elist_JLL[:,1] ]  .- [Z[j]  for j in elist_JLL[:,2] ] ) .* ( [Z[j]  for j in elist_JLL[:,3] ]  .- [Z[j]  for j in elist_JLL[:,4] ] )
            Bii_second_movers = Bii_second_movers .+ ( [ZB[j]  for j in elist_JLL[:,1] ]  .- [ZB[j]  for j in elist_JLL[:,2] ] ) .* ( [ZB[j]  for j in elist_JLL[:,3] ]  .- [ZB[j]  for j in elist_JLL[:,4] ] )

            if settings.first_id_effects == true
                Bii_first_movers = Bii_first_movers .+  ( [ZB_first[j]  for j in elist_JLL[:,1] ]  .- [ZB_first[j]  for j in elist_JLL[:,2] ] ) .* ( [ZB_first[j]  for j in elist_JLL[:,3] ]  .- [ZB_first[j]  for j in elist_JLL[:,4] ] )
            end

            if settings.cov_effects == true
                ZB[N+1:end] .= ZB[N+1:end]
                ZB_first[N+1:end] .= ZB_first[N+1:end]
                Bii_cov_movers = Bii_cov_movers .+ ( [ZB[j]  for j in elist_JLL[:,1] ]  .- [ZB[j]  for j in elist_JLL[:,2] ] ) .* ( [ZB_first[j]  for j in elist_JLL[:,3] ]  .- [ZB_first[j]  for j in elist_JLL[:,4] ] )
                # Bii_cov_movers = Bii_cov_movers .+ ( [ZB[j]  for j in elist_JLL[:,1] ]  .- [ZB[j]  for j in elist_JLL[:,2] ] ) .* ( [ZB_first[j]  for j in elist_JLL[:,3] ]  .- [ZB_first[j]  for j in elist_JLL[:,4] ] )
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

    #Lambda B cov(fe,pe)
    if settings.cov_effects == true
        Lambda_B_cov=sparse(rows,cols,Bii_cov,NT,NT)
        Lambda_B_cov=Lambda_B_cov+triu(Lambda_B_cov,1)'
    end

    #Lambda B, var(pe)
    if settings.first_id_effects == true
        Lambda_B_first=sparse(rows,cols,Bii_first,NT,NT)
        Lambda_B_first=Lambda_B_first+triu(Lambda_B_first,1)'
    end


    if settings.first_id_effects == false & settings.cov_effects == false
        return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second)
    elseif settings.first_id_effects == true & settings.cov_effects == false
        return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first)
    elseif settings.first_id_effects == true  & settings.cov_effects == true
        return (Lambda_P = Lambda_P, Lambda_B_second=Lambda_B_second, Lambda_B_first=Lambda_B_first, Lambda_B_cov=Lambda_B_cov)
    end

end



#9) Creates match first_id using second_id id
function compute_matchid(second_id,first_id)
    match_id2 = string.(first_id).*string.("+",second_id)
    match_id2 = indexin(match_id2, unique(match_id2))

    return match_id2
end

"""
    leave_out_estimation(y,first_id,second_id,controls,settings)

Returns the bias-corrected components, the vector of coefficients, the corresponding fixed effects for every observation, and the diagonal matrices containing the Pii and Biis. At the current version only `controls=nothing' is supported.
"""
#10) Leave Out Function
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
    println("Plug-in Variance of second Effects: ", σ2_ψ_AKM )
    σ2_α_AKM = var(pe)
    println("Plug-in Variance of first Effects: ", σ2_α_AKM )
    σ2_ψα_AKM = cov(pe,-fe)
    println("Plug-in Covariance of second-first Effects: ", σ2_ψα_AKM, "\n")


    #Part 2: Compute Pii, Bii
    @unpack Lambda_P, Lambda_B_second,Lambda_B_first, Lambda_B_cov = eff_res(settings.leverage_algorithm, X,first_id,second_id,match, K, settings)

    #Compute Leaveo-out residual
    I_Lambda_P = I-Lambda_P
    eta_h = I_Lambda_P\eta

    #Compute bias corrected variance comp of second (Firm) Effects
    θ_second = σ2_ψ_AKM -(1/NT)*y'*Lambda_B_second*eta_h

    θ_first = settings.first_id_effects==true  ? σ2_α_AKM -(1/NT)*y'*Lambda_B_first*eta_h : nothing

    θCOV = settings.cov_effects==true ? σ2_ψα_AKM -(1/NT)*y'*Lambda_B_cov*eta_h : nothing


    return (θ_first = θ_first, θ_second = θ_second, θCOV = θCOV, β = beta, Dalpha = pe, Fpsi = fe, Pii = diag(Lambda_P), Bii_first = diag(Lambda_B_first),
            Bii_second = diag(Lambda_B_second), Bii_cov = diag(Lambda_B_cov))

end


"""
    get_leave_one_out_set(y, first_id, second_id, settings, controls)

Returns a tuple with the observation number of the original dataset that belongs to the Leave-out connected set as described in Kline,Saggio, Solvesten. It also provides the corresponding outcome and identifiers in this connected set. At the current version only `controls=nothing' is supported.
"""
function get_leave_one_out_set(y, first_id, second_id, settings, controls)
    @assert settings.first_id_effects == true && settings.cov_effects == true

    # compute y, id firmid, controls, settings
    # compute y, first_id second_id, controls, settings
    (settings.print_level > 0) && println("Finding the leave-one-out connected set")
    @unpack obs_id,  y  , first_id , second_id  = find_connected_set(y,first_id,second_id,settings)
    @unpack obs_id,  y  , first_id , second_id  = prunning_connected_set(y,first_id,second_id, obs_id,settings)
    @unpack obs_id,  y  , first_id , second_id  = drop_single_obs(y,first_id,second_id, obs_id)
    controls == nothing ? nothing : controls = controls[obs_id,:]

    return (obs = obs_id, y = y, first_id = first_id, second_id = second_id, controls = controls)
end

# Do everything naively with no inplace operations, just to get the desired result
# function compute_whole(y,first_id,second_id,controls,settings::VCHDFESettings)

#     @assert settings.first_id_effects == true && settings.cov_effects == true

#     # compute y, id firmid, controls, settings
#     # compute y, first_id second_id, controls, settings
#     (settings.print_level > 0) && println("Finding the leave-one-out connected set")
#     @unpack obs_id,  y  , first_id , second_id  = find_connected_set(y,first_id,second_id;settings)
#     @unpack obs_id,  y  , first_id , second_id  = prunning_connected_set(y,first_id,second_id, obs_id;settings)
#     @unpack obs_id,  y  , first_id , second_id  = drop_single_obs(y,first_id,second_id, obs_id)
#     controls == nothing ? nothing : controls = controls[obs_id,:]

#     if settings.print_level > 0
#         # compute the number of movers
#         num_movers = length(unique(compute_movers(first_id,second_id).movers .* first_id)) - 1 

#         println("\n","Summary statistics of the leave-one-out connected set:")
#         println("Number of observations: ", length(obs_id))
#         println("Number of $(settings.first_id_display_small)s: ", length(unique(first_id)))
#         println("Number of movers: ", num_movers)
#         println("Mean $(settings.observation_id_display_small): ", mean(y))
#         println("Variance of $(settings.observation_id_display_small): ", var(y))
#     end

#     @unpack θ_first, θ_second, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,first_id,second_id,controls,settings)

#     if settings.print_level > 0
#         println("Bias-Corrected Variance Components:")
#         println("Bias-Corrected variance of $(settings.first_id_display_small): $θ_first")
#         println("Bias-Corrected variance of $(settings.second_id_display_small): $θ_second")
#         println("Bias-Corrected covariance of $(settings.first_id_display_small)-$(settings.second_id_display_small) effects: $θCOV")
#     end

#     return (θ_first = θ_first, θ_second = θ_second, θCOV = θCOV, obs = obs_id, β = β, Dalpha = Dalpha, Fpsi = Fpsi, Pii = Pii, Bii_first = Bii_first,
#             Bii_second = Bii_second, Bii_cov = Bii_cov)
# end

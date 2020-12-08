#5) Finds observations connected for every cluster

function check_clustering(clustering_var)
    NT = length(clustering_var)
    counter = ones(size(clustering_var,1));

    #Indexing obs number per worked/id
    gcs = Int.(@transform(groupby(DataFrame(counter = counter, clustering_var = clustering_var), :clustering_var), gcs = cumsum(:counter)).gcs)
    maxD = maximum(gcs)

    index = collect(1:NT)

    #This will be the output, and we will append observations to it in the loop
    list_final=DataFrame(row = Int64[],col = Int64[], id_cluster = Int64[])

    for t=1:maxD
        rowsel =  findall(x->x==true,gcs.==t)
        list_base = DataFrame( id_cluster= [clustering_var[x] for x in rowsel], row =
        [index[x] for x in rowsel]  )

            for tt=t:maxD
                colsel =  findall(x->x==true,gcs.==tt)
                list_sel =  DataFrame( id_cluster= [clustering_var[x] for x in colsel], col =
                [index[x] for x in colsel] )

                merge = outerjoin(list_base, list_sel, on = :id_cluster)
                merge = dropmissing(merge)
                merge = merge[:,[:row, :col, :id_cluster]]
                append!(list_final, merge)

            end
        end

        sort!(list_final, (:row))
        list_final = Matrix(list_final)

        return (list_final = list_final , nnz_2 = size(list_final,1) )

end

# Assumes that obs_leaveoutset is an array of zeros
function whole_prunning!(obs_leaveoutset,y,first_id,second_id;settings)
    obs,  y  , first_id , second_id  = find_connected_set(y,first_id,second_id,settings)
    obs,  y  , first_id , second_id  = prunning_connected_set(y,first_id,second_id,obs,settin)
    obs,  y  , first_id , second_id  = drop_single_obs(y,first_id,second_id,obs)
    for x in obs
        obs_leaveoutset[x] = 1
    end
    return nothing
end



#7) Computes exact Pii to do inference (may be deprecated in the future)
function do_Pii(X, clustering_var)
    n=size(X,1)

    #If no clustering, Lambda_P is just diagonal matrix.
    if clustering_var == nothing
        clustering_var = collect(1:n)
    end

    #Set matrices for parallel environment.
    xx=X'*X
    #P = aspreconditioner(ruge_stuben(xx))
    ldli, la = computeLDLinv(xx)
    buffs = zeros(size(la)[1],Threads.nthreads())
    compute_sol = []
    for i in 1:Threads.nthreads()
        P = approxcholOperator(ldli,buffs[:,i])
        push!(compute_sol,approxcholSolver(P,la))
    end

    #Return the structure of the indexes associated with the clustering variable
    elist, nnz_2 = check_clustering(clustering_var)
    M=size(elist,1)

    #Set elist
    elist_1 = elist[:,1]
    elist_2 = elist[:,2]
    Pii=zeros(M)

    Threads.@threads for i=1:M
        #zexact = zeros(size(X,2))
        col = elist_2[i]
        row = elist_1[i]
        #cg!(zexact, xx, X[col,:], Pl = P , log=true, maxiter=300)
        zexact = compute_sol[Threads.threadid()]([X[col,:]...],verbose=false)

        Pii[i]= (SparseMatrixCSC{Float64,Int64}(X[col,:])'*zexact)[1]
    end

    Lambda_P = sparse(elist[:,1],elist[:,2],Pii,n,n)

    return Lambda_P
end


#8) Performs Statistical Inference on Results
function lincom_KSS(y,X,Z,Transform,clustering_var,Lambda_P; joint_test =false, labels=nothing, restrict=nothing, nsim = 10000, settings = settings)
    #SET DIMENSIONS
    n=size(X,1)
    K=size(X,2)
    #Add Constant
    Z=hcat(ones(size(Z,1)), Z)

    # PART 1: ESTIMATE HIGH DIMENSIONAL MODEL
    xx=X'*X
    xy=X'*y
    compute_sol = approxcholSolver(xx;verbose = settings.print_level > 0)
    beta = compute_sol([xy...];verbose=false)
    eta=y-X*beta

    # PART 1B: VERIFY LEAVE OUT COMPUTATION
    if Lambda_P == nothing
        Lambda_P=do_Pii(X,clustering_var)
    end

    if Lambda_P != nothing && clustering_var !=nothing
        nnz_1=nnz(Lambda_P)
        nnz_2=check_clustering(clustering_var).nnz_2

        if nnz_1 == nnz_2
            println("The structure of the specified Lambda_P is consistent with the level of clustering required by the user.")
        elseif nnz_1 != nnz_2
            error("The user wants cluster robust inference but the Lambda_P provided by the user is not consistent with the level of clustering asked by the user. Try to omit input Lambda_P when running lincom_KSS")
        end
    end
    I_Lambda_P = I-Lambda_P
    eta_h = I_Lambda_P\eta

    #PART 2: SET UP MATRIX FOR SANDWICH FORMULA
    rows,columns, V = findnz(Lambda_P)

    aux= 0.5*(y[rows].*eta_h[columns] + y[columns].*eta_h[rows])
    sigma_i=sparse(rows,columns,aux,n,n)

    aux= 0.5*(eta[rows].*eta[columns] + eta[columns].*eta[rows])
    sigma_i_res=sparse(rows,columns,aux,n,n)

    r=size(Z,2);
    wy=Transform*beta
    zz=Z'*Z

    numerator=Z\wy
    chet=wy-Z*numerator
    aux= 0.5*(chet[rows].*chet[columns] + chet[columns].*chet[rows])
    sigma_i_chet=sparse(rows,columns,aux,n,n)

    #PART 3: COMPUTE
    denominator=zeros(r,1)
    denominator_RES=zeros(r,1)

    for q=1:r
        v=sparse([q],[1.0],[1.0],r,1)
        v=zz\[v...]
        v=Z*v
        v=Transform'*v

        right = compute_sol(v;verbose=false)

        left=right'

        denominator[q]=left*(X'*sigma_i*X)*right
        denominator_RES[q]=left*(X'*sigma_i_res*X)*right
    end

    test_statistic=numerator./(sqrt.(denominator))
    #zz_inv=zz^(-1)
    SE_linear_combination_NAI=zz\(Z'*sigma_i_chet*Z)/zz


    #PART 4: REPORT
    println("Inference on Linear Combinations:")
    if labels == nothing
        for q=2:r
            if q <= r
                println("Linear Combination - Column Number ", q-1," of Z: ", numerator[q] )
                println("Standard Error of the Linear Combination - Column Number ", q-1," of Z: ", sqrt(denominator[q]) )
                println("T-statistic - Column Number ", q-1, " of Z: ", test_statistic[q])
            end
        end
    else
        for q=2:r
            tell_me = labels[q-1]
            println("Linear Combination associated with ", tell_me,": ", numerator[q] )
            println("Standard Error  associated with ", tell_me,": ", sqrt(denominator[q]) )
            println("T-statistic  associated with ", tell_me,": ", test_statistic[q])
        end
    end


    # PART 5: Joint-test. Quadratic form beta'*A*beta
    if joint_test == true

        if restrict  == nothing
            restrict=sparse(collect(1:r-1),collect(2:r),1.0,r-1,r)
        end

        v=restrict*(zz\(Z'*Transform))
        v=v'
        #v=sparse(v) #ldiv doesn't work for sparse RHS
        r=size(v,2)

        #Auxiliary
        aux=xx\v[:,:]
        opt_weight=v'*aux
        opt_weight=opt_weight^(-1)
        opt_weight=(1/r)*(opt_weight+opt_weight')/2

        #Eigenvalues, eigenvectors, and relevant components
        lambda , Qtilde = eigs( v*opt_weight*v', xx; nev=r,ritzvec=true)
        lambda = Real.(lambda)
        Qtilde = Real.(Qtilde)
        #lambdaS, QtildeS = eigs(v*opt_weight*v', xx; nev=1,which=:SM,ritzvec=true)
        #lambdaS = [lambdaL; lambdaS]
        #Qtilde = hcat(QtildeL,QtildeS)

        W=X*Qtilde
        V_b=W'*sigma_i*W
        V_b = (1/2)*(V_b + V_b')

        #Now focus on obtaining matrix Lambda_B with the A test associated with a joint hypothesis testing.
        Bii=opt_weight^(0.5)*aux';
        Bii=Bii*X'
        Bii=Bii'
        Bii = 0.5*(Bii[rows,:].*Bii[columns,:] + Bii[columns,:].*Bii[rows,:])
        Bii = sum(Bii,dims=2)[:]
        Lambda_B=sparse(rows,columns,Bii,n,n)

        #Leave Out Joint-Statistic
        stat=(v'*beta)'*opt_weight*(v'*beta)-y'*Lambda_B*eta_h

        #Now simulate critical values under the null.
        mu=zeros(r)
        sigma = V_b
        b_sim = MvNormal(mu,sigma)
        b_sim = rand(b_sim, nsim)

        #theta_star_sim=sum(lambda'.*(b_sim.^2 - diag(V_b)'),2)
        theta_star_sim = sum(lambda.*b_sim.^2 .- lambda.* diag(V_b),dims=1)
        pvalue=mean(theta_star_sim.>stat)

        #Report
        println("Joint-Test Statistic: ", stat)
        println("p-value: ", pvalue)

    end

    test_statistic=test_statistic[2:end]
    linear_combination=numerator[2:end]
    SE_linear_combination_KSS=sqrt.(denominator[2:end])
    SE_linear_combination_RES=sqrt.(denominator_RES[2:end])
    SE_linear_combination_NAI=diag(SE_linear_combination_NAI)
    SE_linear_combination_NAI=sqrt.(SE_linear_combination_NAI[2:end])

    return nothing
end

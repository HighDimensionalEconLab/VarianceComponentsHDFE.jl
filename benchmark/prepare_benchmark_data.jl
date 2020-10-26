using SparseArrays, LinearAlgebra
using JLD, CSV, DataFrames, DataFramesMeta
#using DataDeps
using ArtifactUtils, Pkg.Artifacts

artifacts_path = joinpath(pkgdir(VarianceComponentsHDFE), "deps\\Artifacts.toml")

# add in code to generate .jld, .csv, etc. in the benchmark/data directory
# This should only generate if the file doesn't exist, or if this `force_generate = true`
force_generate = false

pkg_dir = pkgdir(VarianceComponentsHDFE)

# Obtain environment variables
run_large_benchmark = get(ENV, "VCHDFE_LARGE_BENCHMARK", "false")  == "true" ? true : false
run_huge_benchmark = get(ENV, "VCHDFE_HUGE_BENCHMARK", "false")  == "true" ? true : false

## Medium-Sized Network Generator
function rademacher!(R; demean = false)
    R .= R - (R .== 0)

    if demean == true
        R .= R .- mean(R, dims = 2)
    end
    return nothing
end

function compute_X_No_Controls(data)
    id = data.id
    firmid = data.firmid
    y = data.y

    NT = size(y,1);
    J = maximum(firmid);
    N = maximum(id);
    K = 0
    nparameters = N + J + K

    #Worker Dummies
    D = sparse(collect(1:NT),id,1);

    #Firm Dummies
    F = sparse(collect(1:NT),firmid,1);

    # N+J x N+J-1 restriction matrix
    S= sparse(1.0I, J-1, J-1);
    S=vcat(S,sparse(-zeros(1,J-1)));
    X_Laplacian = hcat(D, -F)
    X_GroundedLaplacian = hcat(D, -F*S)
    S_xx = Symmetric(X_GroundedLaplacian'*X_GroundedLaplacian)
    S_xx_lap = Symmetric(SparseMatrixCSC{Float64,Int64}(X_Laplacian'*X_Laplacian))

    A, diag = adj(SparseMatrixCSC{Float64,Int64}(X_Laplacian'*X_Laplacian))

    return X_Laplacian, X_GroundedLaplacian, S_xx, S_xx_lap, A
end

function compute_X_Controls(data)
    id = data.id
    firmid = data.firmid
    y = data.y


    NT=size(y,1);
    J=maximum(firmid);
    N=maximum(id);
    K = 2
    nparameters = N + J + K

    #Worker Dummies
    D=sparse(collect(1:NT),id,1);

    #Firm Dummies
    F=sparse(collect(1:NT),firmid,1);

    # N+J x N+J-1 restriction matrix
    S= sparse(1.0I, J-1, J-1);
    S=vcat(S,sparse(-zeros(1,J-1)));

    #Assuming columns 5 and 6 are the controls
    controls = hcat(data.control1[id], data.control2[id])

    Xcontrols = hcat(D,F*S,controls)
    S_xx = Symmetric(Xcontrols'*Xcontrols)

    return Xcontrols, S_xx
end

# We don't use the symmetric flag here, it is making the cluster kill the process
# Also flagged X_laplacian as a sparse matrix of Floats
function compute_X_No_Controls_Huge(data)
    id = data.id
    firmid = data.firmid
    y = data.y

    NT = size(y,1);
    J = maximum(firmid);
    N = maximum(id);
    K = 0
    nparameters = N + J + K

    #Worker Dummies
    D = sparse(collect(1:NT),id,1);

    #Firm Dummies
    F = sparse(collect(1:NT),firmid,1);

    # N+J x N+J-1 restriction matrix
    S= sparse(1.0I, J-1, J-1);
    S=vcat(S,sparse(-zeros(1,J-1)));
    X_Laplacian = SparseMatrixCSC{Float64,Int64}(hcat(D, -F))
    X_GroundedLaplacian = hcat(D, -F*S)
    S_xx = X_GroundedLaplacian'*X_GroundedLaplacian
    S_xx_lap = X_Laplacian'*X_Laplacian
    A, diag = adj(S_xx_lap)

    return X_Laplacian, X_GroundedLaplacian, S_xx, S_xx_lap, A
end

if ~isfile(pkg_dir*"/benchmark/data/medium_main.jld") || force_generate
    #data = CSV.read(datadep"VarianceComponentsHDFE/medium_nocontrols_pruned.csv"; header=true)
    data_hash = artifact_hash("medium_nocontrols_pruned",artifacts_path)
    data = CSV.read(joinpath(artifact_path(data_hash), "medium_nocontrols_pruned.csv"))

    X_Laplacian, X_GroundedLaplacian, S_xx , S_xx_lap, A = compute_X_No_Controls(data)
    save(pkg_dir*"/benchmark/data/medium_main.jld", "X_Laplacian", X_Laplacian, "X_GroundedLaplacian", X_GroundedLaplacian, "S_xx", S_xx, "S_xx_lap", S_xx_lap, "A", A)
end

if ~isfile(pkg_dir*"/benchmark/data/medium_controls_main.jld") || force_generate
    #data = CSV.read(datadep"VarianceComponentsHDFE/medium_controls_pruned.csv"; header=true)
    data_hash = artifact_hash("medium_controls_pruned",artifacts_path)
    data = CSV.read(joinpath(artifact_path(data_hash), "medium_controls_pruned.csv"))

    Xcontrols, S_xx = compute_X_Controls(data)
    save(pkg_dir*"/benchmark/data/medium_controls_main.jld", "Xcontrols", Xcontrols, "S_xx", S_xx)
end

if run_large_benchmark && (~isfile(pkg_dir*"/benchmark/data/large_main.jld") || force_generate)
    #data = CSV.read(datadep"VarianceComponentsHDFE/large_nocontrols_pruned.csv"; header=true)
    data_hash = artifact_hash("large_nocontrols_pruned",artifacts_path)
    data = CSV.read(joinpath(artifact_path(data_hash), "large_nocontrols_pruned.csv"))

    X_Laplacian, X_GroundedLaplacian, S_xx , S_xx_lap, A = compute_X_No_Controls(data)
    save(pkg_dir*"/benchmark/data/large_main.jld", "X_Laplacian", X_Laplacian, "X_GroundedLaplacian", X_GroundedLaplacian, "S_xx", S_xx, "S_xx_lap", S_xx_lap, "A", A)
 end

 if run_large_benchmark && (~isfile(pkg_dir*"/benchmark/data/large_controls_main.jld") || force_generate)
    #lar_csv = CSV.read(joinpath(artifact"large_nocontrols_pruned", "large_nocontrols_pruned.csv"))
    data_hash = artifact_hash("large_controls_pruned",artifacts_path)
    data = CSV.read(joinpath(artifact_path(data_hash), "large_controls_pruned.csv"))

    Xcontrols, S_xx = compute_X_Controls(data)
    save(pkg_dir*"/benchmark/data/large_controls_main.jld", "Xcontrols", Xcontrols, "S_xx", S_xx)
 end

 if run_huge_benchmark && (~isfile(pkg_dir*"/benchmark/data/huge_main.jld") || force_generate)
     #data = CSV.read(datadep"VarianceComponentsHDFE/huge_pruned_main.csv"; header=true)
     data_hash = artifact_hash("huge_nocontrols_pruned",artifacts_path)
     data = CSV.read(joinpath(artifact_path(data_hash), "huge_nocontrols_pruned.csv"))

     X_Laplacian, X_GroundedLaplacian, S_xx , S_xx_lap, A = compute_X_No_Controls_Huge(data)
     save(pkg_dir*"/benchmark/data/huge_main.jld", "X_Laplacian", X_Laplacian, "X_GroundedLaplacian", X_GroundedLaplacian, "S_xx", S_xx, "S_xx_lap", S_xx_lap, "A", A)
 end
 
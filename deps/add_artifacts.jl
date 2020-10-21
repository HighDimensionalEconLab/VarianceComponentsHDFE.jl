#  NOTE:  This only needs to be run when the binaries change or need to be modified.
using ArtifactUtils, Pkg.Artifacts, VarianceComponentsHDFE, CSV
artifacts_path = joinpath(pkgdir(VarianceComponentsHDFE), "Artifacts.toml")
url_prefix = "https://vchdfe.s3-us-west-2.amazonaws.com"
add_artifact!(
            artifacts_path,
           "medium_nocontrols_pruned",
           "$url_prefix/medium_pruned_main.tar.gz",
           force=true,
           lazy = true
       )

add_artifact!(
    artifacts_path,
    "medium_controls_pruned",
    "$url_prefix/medium_controls_pruned_main.tar.gz",
    force=true,
    lazy = true
)


add_artifact!(
    artifacts_path,
    "large_nocontrols_pruned",
    "$url_prefix/large_pruned_main.tar.gz",
    force=true,
    lazy =true
)

add_artifact!(
    artifacts_path,
    "large_controls_pruned",
    "$url_prefix/large_controls_pruned_main.tar.gz",
    force=true,
    lazy = true
)

add_artifact!(
    artifacts_path,
    "huge_nocontrols_pruned",
    "$url_prefix/huge_pruned_main.tar.gz",
    force=true,
    lazy = true
)

#Ensure Installation of the files 

ensure_artifact_installed("medium_nocontrols_pruned",artifacts_path)
ensure_artifact_installed("medium_controls_pruned",artifacts_path)
ensure_artifact_installed("large_nocontrols_pruned",artifacts_path)
ensure_artifact_installed("large_controls_pruned",artifacts_path)
ensure_artifact_installed("huge_nocontrols_pruned",artifacts_path)

# Testing 
med_csv_hash = artifact_hash("medium_nocontrols_pruned",artifacts_path)
med_csv = CSV.read(joinpath(artifact_path(med_csv_hash), "medium_nocontrols_pruned.csv"))

lar_csv_hash = artifact_hash("large_nocontrols_pruned",artifacts_path)
lar_csv = CSV.read(joinpath(artifact_path(lar_csv_hash), "large_nocontrols_pruned.csv"))

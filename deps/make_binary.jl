using PackageCompiler
create_app(".", "./vchdfe", app_name="vchdfe", force=true, precompile_execution_file="deps/precompilation_execution_file.jl", include_lazy_artifacts = false)

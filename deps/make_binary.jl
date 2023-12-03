using PackageCompiler
create_app(".", "./vchdfe", executables=["vchdfe"=>"julia_main"], force=true, precompile_execution_file="deps/precompilation_execution_file.jl", include_lazy_artifacts=true)

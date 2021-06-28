using CSV
using VarianceComponentsHDFE
using Parameters
using SparseArrays
using LinearAlgebra
using DataFrames
using Pipe
using ShiftedArrays
using LinearMaps
using Random
using Statistics

# data = CSV.read("data/gen_data_homoskedastic_v03_v01_psibar.csv", DataFrame, missingstrings = ["NA", ""])
# data = CSV.read("data/DGP1_no_error2.csv", DataFrame, missingstrings = ["NA", ""])
data = CSV.read("data/first_step_reduced_leave_out_data_renamed.csv", DataFrame, missingstrings = ["NA", ""])




function weightedmean(x, w)
    wxsum = wsum = 0
    for (x,w) in zip(x,w)
        wx = w*x
        if !ismissing(wx)
            wxsum += wx
            wsum += w
        end
    end
    return wxsum / wsum
end
## Binned Scatter plot
# ] add Binscatters
using RDatasets, Plots, StatsModels, CategoricalArrays
df_bs = @pipe groupby(data, [:firmidg, :year]) |> combine(_, :psi => first => :psi, :alpha => mean => :alphabar, nrow => :nworkers)
sort(data, [:firmidg, :year])
df_bs2 = df_bs
df_bs2.:psi_d = df_bs2.:psi
df_bs2.:alpha_d = df_bs2.alphabar
# df_bs2 = @pipe groupby(df_bs, [:firmidg]) |> transform(_, :psi => (x -> (x .- lag(x)) ) => :psi_d, :alphabar => (x -> (x .- lag(x))) => :alpha_d)
df_bs2 = dropmissing(df_bs2)
# binscatter(df_bs2, @formula(psi_d ~ alpha_d), weights = :nworkers, 100, seriestype = :linearfit; title = "binned scatter plot", label = "PI")


x = df_bs2.alpha_d
y = df_bs2.psi_d
w = df_bs2.nworkers

x = vcat(fill.(x, w)...)
y = vcat(fill.(y, w)...)

df = DataFrame(x = x, y = y)

xtile(x; n=4) = levelcode.(cut(x, n; allowempty = true))

df = @pipe transform(df, :x => (x -> xtile(x; n=90)) => :q)

df = @pipe sort(df, :x) |> groupby(_, :q) |> combine(_, :x => mean, :y => mean)


const_a_kss = weightedmean(df_bs2.:psi_d - (-0.05614690610628232 * df_bs2.:alpha_d), df_bs2.:nworkers)

f_kss(x) = const_a_kss + (-0.05614690610628232 * x)

const_a_PI = weightedmean(df_bs2.psi_d - (-0.08072056306897829 * df_bs2.:alpha_d), df_bs2.:nworkers)
f_PI(x) = const_a_PI + (-0.08072056306897829 * x)
lowend = -0.31
highend = 0.22

plot(df.:x_mean, df.:y_mean, seriestype = :scatter;  title = "Binned scatter plot", label = "mean")

plot!(f_kss, lowend, highend; label = "KSS", legend = true)
plot!(f_PI, lowend, highend; label = "PI")
plot!(size = (900, 680))
xlabel!("change in person effect")
ylabel!("change in firm-year effect")
fontsize = 12
kss_string = "y = 0.017 - 0.056x \n R2 = 0.004 "
annotate!(-0.18, 0.026, text(kss_string, :right, fontsize, color = :brown))
pi_string = "y = 0.017 - 0.081x \n R2 = 0.008 \n SE = 0.016 "
annotate!(-0.13, 0.036, text(pi_string, :right, fontsize, color = :green))


png("binscatterplot2")

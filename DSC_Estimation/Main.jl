using JLD2
using CSV
using Random
using Dates
using LinearAlgebra
using Statistics

# Data Structure
include("InsChoiceData.jl")

#Halton Draws
include("Halton.jl")

# Random Coefficients MLE
include("RandomCoefficients.jl")
include("RandomCoefficients_der.jl")
include("Log_Likehood.jl")
include("Estimate_Basic.jl")
include("utility.jl")
include("Specification_Run.jl")
println("Code Loaded")

# Load the Data
include("load.jl")
df = df[df[:gra].==16,:]
println("Data Loaded")

#### General Specification ####
spec_per = [:hh_year_id]
spec_prd = [:product]
spec_ch = [:choice]
spec_ch_last = [:iplan]
spec_wgt=[:constant]
mixed_draws = 100

### interact default plans and "active" choice with preference parameters

rundate = Dates.today()

#### Mixed Logit DSC Specifications ####
println("###############################")
println("Mixed Specification 1")
filename = "ML_spec1_$rundate"
out = MainSpec(df,filename,
    haltonDim = 500,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= [:padj,:iplan],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:active,:dprem,:def_padj, :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,:def_mtl_plt,:def_mtl_s73,:def_mtl_s87,:def_mtl_s94],
    spec_demR=Vector{Symbol}(undef,0),
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method="ga",
    ga_itr = 200)


#### Basic Switching Cost Specifications ####
println("###############################")
println("Basic Specification 1")
filename = "basic_spec1_$rundate"
df_est = df[(df[:hasi].==1).&(df[:dcat].==1),:]
basic_out_1 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=Vector{Symbol}(undef,0),
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt)

println("###############################")
println("Basic Specification 2")
filename = "basic_spec2_$rundate"
basic_out_2 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=Vector{Symbol}(undef,0),
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:product],
    spec_wgt= spec_wgt)


println("###############################")
println("Basic Specification 3")
filename = "basic_spec3_$rundate"
basic_out_3 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt)

println("###############################")
println("Basic Specification 4")
filename = "basic_spec4_$rundate"
basic_out_4 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:product],
    spec_wgt= spec_wgt)

#### Basic Searching/Switching Specifications ####

println("###############################")
println("Basic Searching Specification 1")
filename = "basic_search_spec1_$rundate"
# x_start, out = basic_out_1
# x_start = vcat([5],rand(6)/10 .-0.05,x_start)
# df_est = df[(df[:hasi].==1).&(df[:dcat].==1),:]
bs_out_1 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    # spec_inertchr= [:constant,:padj],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:active,:dprem],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method="ga")


println("###############################")
println("Basic Searching Specification 2")
filename = "basic_search_spec2_$rundate"
bs_out_2 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:active,:dprem],
    spec_demR=Vector{Symbol}(undef,0),
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:product],
    spec_wgt= spec_wgt,
    method="ga")


println("###############################")
println("Basic Searching Specification 3")
filename = "basic_search_spec3_$rundate"
# x_start, out = basic_out_3
# x_start = vcat([20],rand(6)/10 .-0.05],x_start)
bs_out_3 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:active,:dprem],
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method = "ga")

println("###############################")
println("Basic Searching Specification 4")
filename = "basic_search_spec4_$rundate"
# x_start, out = basic_out_4
# x_start = vcat([20],rand(6)/10 .-0.05],x_start)
bs_out_4 = MainSpec(df,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr=[:dprem,:active],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:active,:dprem],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:product],
    spec_wgt= spec_wgt,
    method = "ga")

#### Mixed Logit Switching Cost Specifications ####
println("###############################")
println("Mixed Specification 1")
filename = "ML_spec1_$rundate"
out = MainSpec(df,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= [:padj,:iplan],
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=Vector{Symbol}(undef,0),
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method="ga",
    ga_itr = 100)

println("###############################")
println("Mixed Specification 2")
filename = "ML_spec2_$rundate"
out = MainSpec(df,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= [:padj,:iplan],
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=Vector{Symbol}(undef,0),
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:product],
    spec_wgt= spec_wgt)


println("###############################")
println("Mixed Specification 3")
filename = "ML_spec3_$rundate"
out = MainSpec(df,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= [:padj,:iplan],
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt)

println("###############################")
println("Mixed Specification 4")
filename = "ML_spec4_$rundate"
out = MainSpec(df,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan],
    spec_prodchr_0= [:padj,:iplan],
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:product],
    spec_wgt= spec_wgt)

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
df_LA = df[df[:gra].==16,:]

### Add Anthem Fixed Effect
df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")

println("Data Loaded")

#### General Specification ####
spec_per = [:hh_year_id]
spec_prd = [:product]
spec_ch = [:choice]
spec_ch_last = [:iplan]
spec_wgt=[:constant]
mixed_draws = 500

### interact default plans and "active" choice with preference parameters

rundate = Dates.today()


# ### Basic Switching Cost Specifications ####
# println("###############################")
# println("Basic Specification 1")
# filename = "basic_spec1_$rundate"
# # df_est = df[(df[:hasi].==1).&(df[:dcat].==1),:]
# basic_out_1 = MainSpec(df,filename,
#     haltonDim = 1,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     spec_inertchr=Vector{Symbol}(undef,0),
#     spec_demR=Vector{Symbol}(undef,0),
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal,:netname],
#     spec_wgt= spec_wgt,
#     method="ga",
#     ga_itr=20)
#
# println("###############################")
# println("Basic Specification 2")
# filename = "basic_spec2_$rundate"
# basic_out_2 = MainSpec(df,filename,
#     haltonDim = 1,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     spec_inertchr=Vector{Symbol}(undef,0),
#     spec_demR=Vector{Symbol}(undef,0),
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:product],
#     spec_wgt= spec_wgt,
#     method="ga",
#     ga_itr=20)
#
#
println("###############################")
println("Basic Specification 3")
filename = "basic_spec3_$rundate"
basic_out_3 = MainSpec(df_LA,filename,
    haltonDim = 1,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,:inet,:iiss],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method="ga",
    ga_itr=20)
#
# println("###############################")
# println("Basic Specification 4")
# filename = "basic_spec4_$rundate"
# basic_out_4 = MainSpec(df,filename,
#     haltonDim = 1,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     spec_inertchr=Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:product],
#     spec_wgt= spec_wgt,
#     method="ga",
#     ga_itr=20)
#
# #### Basic Searching/Switching Specifications ####
#
# println("###############################")
# println("Basic Searching Specification Test")
# filename = "basic_search_test_$rundate"
# # x_start, out = basic_out_1
# # x_start = vcat([5],rand(6)/10 .-0.05,x_start)
# df_est = df[(df[:hasi].==1).&(df[:dcat].==1),:]
# bs_out_1 = MainSpec(df,filename,
#     haltonDim = 1,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     # spec_inertchr= [:constant,:padj],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
#                     :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
#                     :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
#                     :def_mtl_s87,:def_mtl_s94],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal,:netname],
#     spec_wgt= spec_wgt,
#     method="ga")
#
# println("###############################")
# println("Basic Searching Specification 1")
# filename = "basic_search_spec1_$rundate"
# bs_out_1 = MainSpec(df,filename,
#     haltonDim = 1,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     # spec_inertchr= [:constant,:padj],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
#                     :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
#                     :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
#                     :def_mtl_s87,:def_mtl_s94],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal,:netname],
#     spec_wgt= spec_wgt,
#     method="ga")
#
# println("###############################")
# println("Basic Searching Specification 2")
# filename = "basic_search_spec2_$rundate"
# bs_out_2 = MainSpec(df,filename,
#     haltonDim = 1,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr =[:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
#                     :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
#                     :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
#                     :def_mtl_s87,:def_mtl_s94],
#     spec_demR=Vector{Symbol}(undef,0),
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:product],
#     spec_wgt= spec_wgt,
#     method="ga")
#
#
println("###############################")
println("Basic Searching Specification 3")
filename = "basic_search_spec3_$rundate"
# x_start, out = basic_out_3
# x_start = vcat([20],rand(6)/10 .-0.05],x_start)
# bs_out_3 = MainSpec(df_LA,filename,
#     haltonDim = 1,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
#                     :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
#                     :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
#                     :def_mtl_s87,:def_mtl_s94],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal,:netname],
#     spec_wgt= spec_wgt,
#     method = "ga")
#
# println("###############################")
# println("Basic Searching Specification 4")
# filename = "basic_search_spec4_$rundate"
# # x_start, out = basic_out_4
# # x_start = vcat([20],rand(6)/10 .-0.05],x_start)
# bs_out_4 = MainSpec(df,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= Vector{Symbol}(undef,0),
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
#                     :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
#                     :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
#                     :def_mtl_s87,:def_mtl_s94],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:product],
#     spec_wgt= spec_wgt,
#     method = "ga")

# ### Mixed Logit Switching Cost Specifications ####
# println("###############################")
# println("Mixed Specification 1")
# filename = "ML_spec1_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= [:padj,:iiss],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal,:netname],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 10)
#
# println("###############################")
# println("Mixed Specification 2")
# filename = "ML_spec2_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= [:padj,:iiss],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
#                     :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
#                     :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
#                     :def_mtl_s87,:def_mtl_s94],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal,:netname],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 100)

# println("###############################")
# println("Mixed Specification 3")
# filename = "ML_spec3_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss],
#     spec_prodchr_0= [:padj,:iiss],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal,:netname],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 10)
#


# println("###############################")
# println("Mixed Specification 5")
# filename = "ML_spec5_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss,
#     :issfe_1, :issfe_2, :issfe_5, :issfe_6,
#     :issfe_8, :issfe_9, # Leave Out LA Care
#     :netfe_2, :netfe_3, :netfe_4, :netfe_7,
#     :netfe_11, :netfe_12, :netfe_13, :netfe_15],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet,:iiss],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 100)

println("###############################")
println("Mixed Specification 6")
filename = "ML_spec6_$rundate"
mx_out_1 = MainSpec(df_LA,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,:inet,:iiss,
    :issfe_1, :issfe_2, :issfe_5, :issfe_6,
    :issfe_8, :issfe_9, # Leave Out LA Care
    :netfe_2, :netfe_3, :netfe_4, :netfe_7,
    :netfe_11, :netfe_12, :netfe_13, :netfe_15],
    spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
                    :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
                    :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
                    :def_mtl_s87,:def_mtl_s94],
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_prodInt=[:padj,:iplan,:inet,:iiss],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 200)

println("###############################")
println("Mixed Specification 4")
filename = "ML_spec4_$rundate"
mx_out_1 = MainSpec(df_LA,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,:inet,:iiss],
    spec_prodchr_0= [:padj,:iplan],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
                    :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
                    :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
                    :def_mtl_s87,:def_mtl_s94],
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 200)


println("###############################")
println("Mixed Specification 7")
filename = "ML_spec7_$rundate"
mx_out_1 = MainSpec(df_LA,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,:inet,:iiss,
    :issfe_1, :issfe_2, :issfe_5, :issfe_6,
    :issfe_8, :issfe_9, # Leave Out LA Care
    :netfe_2, :netfe_3, :netfe_4, :netfe_7,
    :netfe_11, :netfe_12, :netfe_13, :netfe_15],
    spec_prodchr_0= [:issfe_1],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
                    :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
                    :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
                    :def_mtl_s87,:def_mtl_s94],
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_prodInt=[:padj,:iplan,:inet,:iiss],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 200)


println("###############################")
println("Mixed Specification 8")
filename = "ML_spec8_$rundate"
mx_out_1 = MainSpec(df_LA,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,:inet,:iiss,
    :issfe_1, :issfe_2, :issfe_5, :issfe_6,
    :issfe_8, :issfe_9, # Leave Out LA Care
    :netfe_2, :netfe_3, :netfe_4, :netfe_7,
    :netfe_11, :netfe_12, :netfe_13, :netfe_15],
    spec_prodchr_0= [:padj],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
                    :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
                    :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
                    :def_mtl_s87,:def_mtl_s94],
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_prodInt=[:padj,:iplan,:inet,:iiss],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal,:netname],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 200)

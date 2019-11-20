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
include("DerivFunctions.jl")
include("Log_Likehood.jl")
include("Estimate_Basic.jl")
include("utility.jl")
include("Specification_Run.jl")
println("Code Loaded")

# Load the Data
include("load.jl")
# df_LA = df[df[:gra].==16,:]
# # df_LA = df[df[:hasi].==1.0,:]
df_LA = df
# ### Add Anthem Fixed Effect
# df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")

println("Data Loaded")

#### General Specification ####
spec_per = [:hh_id]
spec_prd = [:product]
spec_ch = [:choice]
spec_ch_last = [:iplan]
spec_wgt=[:constant]
mixed_draws = 500

### interact default plans and "active" choice with preference parameters

rundate = Dates.today()
# file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/ML_spec6_2019-05-25.jld2"
# @load file p_est spec_Dict fval


[:padj,:iplan,
:issfe_1, :issfe_2, :issfe_3, :issfe_4,
:issfe_6, :issfe_7, # Leave Out LA Care
:netfe_2, :netfe_3, :netfe_5, :netfe_6,
:netfe_8, :netfe_9],

println("###############################")
println("Specification 1")
### Only plan-level switching cost, no Inertia
filename = "Spec1_$rundate"
mx_out_1 = MainSpec(df_LA,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,
    :issfe_1, :issfe_2, :issfe_3, :issfe_4,
    :issfe_6, :issfe_7, # Leave Out LA Care
    :netfe_2, :netfe_3, :netfe_5, :netfe_6,
    :netfe_8, :netfe_9],
    spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_3, :issfe_4],
    spec_inertchr= Vector{Symbol}(undef,0),
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_prodInt=[:padj,:iplan],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 50,ll_start = true)

#
println("###############################")
println("Specification 2")
### Fully Specified Switching Costs, no Inertia
filename = "Spec2_$rundate"
mx_out_1 = MainSpec(df_LA,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,
    :issfe_1, :issfe_2, :issfe_3, :issfe_4,
    :issfe_5, :issfe_6, :issfe_8, # Leave Out LA Care
    :issfe_9, :issfe_10, :issfe_11, :issfe_12, :issfe_13,
    :netfe_2, :netfe_3, :netfe_5, :netfe_6,
    :netfe_10, :netfe_11, :netfe_12],
    spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
    spec_inertchr= Vector{Symbol}(undef,0),
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_prodInt=[:padj,:iplan,:inet],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 50,ll_start = true)

println("###############################")
println("Specification 3")
### Full Specification
filename = "Spec3_$rundate"
mx_out_1 = MainSpec(df_LA,filename,
    haltonDim = mixed_draws,
    spec_per = spec_per,
    spec_prd = spec_prd,
    spec_ch = spec_ch,
    spec_ch_last = spec_ch_last,
    spec_prodchr = [:padj,:iplan,
    :issfe_1, :issfe_2, :issfe_3, :issfe_4,
    :issfe_5, :issfe_6, :issfe_8, # Leave Out LA Care
    :issfe_9, :issfe_10, :issfe_11, :issfe_12, :issfe_13,
    :netfe_2, :netfe_3, :netfe_5, :netfe_6,
    :netfe_10, :netfe_11, :netfe_12],
    spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:autodp,
                        #Metal Fixed Effects
                        :def_mtl_cat,:def_mtl_gld, # Leave Out Bronze
                        :def_mtl_hdp,:def_mtl_plt,
                        :def_mtl_slv,:def_mtl_s73,
                        :def_mtl_s87,:def_mtl_s94,
                        # Network Fixed Effects
                        # :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
                        # :def_issfe_6, :def_issfe_7, # Leave Out LA Care
                        # :def_netfe_4, :def_netfe_6, # Drop net02, net03
                        # :def_netfe_9, :def_netfe_12, # Drop net10, net08
                        # Year Fixed Effects
                        :year_2016,:year_2017,:year_2018,
                        #Rating Area Fixed Effect
                        ],
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_prodInt=[:padj,:iplan,:inet],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 50,ll_start=true)


#
# println("###############################")
# println("Specification 2, with Kaiser Robust Check")
# ### Fully Specified Switching Costs, no Inertia
# filename = "Spec2_kais_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss,:iplan_kais,:inet_kais,
#     :issfe_1, :issfe_2, :issfe_3, :issfe_4,
#     :issfe_6, :issfe_7, # Leave Out LA Care
#     :netfe_2, :netfe_3, :netfe_4, :netfe_6,
#     :netfe_8, :netfe_9, :netfe_10, :netfe_12],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_3, :issfe_4],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet,:iiss,:iplan_kais,:inet_kais],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)
#
#
# println("###############################")
# println("Specification 3, with Kaiser Robust Check")
# ### Full Specification
# filename = "Spec3_kais_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:iiss,:iplan_kais,:inet_kais,
#     :issfe_1, :issfe_2, :issfe_3, :issfe_4,
#     :issfe_6, :issfe_7, # Leave Out LA Care
#     :netfe_2, :netfe_3, :netfe_4, :netfe_6,
#     :netfe_8, :netfe_9, :netfe_10, :netfe_12],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_3, :issfe_4],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,
#                         #Metal Fixed Effects
#                         :def_mtl_brz,:def_mtl_cat,:def_mtl_gld, # Leave Out Silver
#                         :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
#                         :def_mtl_s87,:def_mtl_s94,
#                         # Network Fixed Effects
#
#                         # :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
#                         # :def_issfe_6, :def_issfe_7, # Leave Out LA Care
#                         # :def_netfe_4, :def_netfe_6, # Drop net02, net03
#                         # :def_netfe_9, :def_netfe_12, # Drop net10, net08
#                         # Year Fixed Effects
#                         :year_2016,:year_2017,:year_2018],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet,:iiss,:iplan_kais,:inet_kais],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)

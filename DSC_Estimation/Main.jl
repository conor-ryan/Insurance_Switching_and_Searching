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

home = "$(homedir())/Documents/Research/CovCAInertia"
# home = "G:/Shared drives/CovCAInertia"

# Load the Data
include("load.jl")
# df_LA = df[df[:gra].==16,:]
# df_LA_active = df_active[df_active[:gra].==16,:]
# # df_LA = df[df[:hasi].==1.0,:]
df_LA = df
df_LA_active = df_active
# ### Add Anthem Fixed Effect
# df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")

println("Data Loaded")

#### General Specification ####
spec_per = [:hh_id]
spec_prd = [:product]
spec_ch = [:choice]
spec_ch_last = [:autoplan]
spec_wgt=[:constant]
mixed_draws = 2000

### interact default plans and "active" choice with preference parameters

rundate = Dates.today()
# file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/ML_spec6_2019-05-25.jld2"
# @load file p_est spec_Dict fval
#
# println("###############################")
# println("Specification 1")
# ### Only plan-level switching cost, no Inertia
# filename = "Spec1_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     # spec_fixEff=[:metal,:netname],
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start = true)
#
#
# println("###############################")
# println("Specification 2")
# ### Fully Specified Switching Costs, no Inertia
# filename = "Spec2_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start = true)
#
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
    spec_prodchr = [:padj,:iplan,:inet],
    spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
    spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,
                        #Continuous Premium
                        # :autodp,
                        # Categorical Premium
                        :autodp_disc_1,:autodp_disc_2,:autodp_disc_4,:autodp_disc_5,:autodp_disc_6,
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
                        #Rating Area Fixed Effect
                        :grafe_2,:grafe_3,:grafe_4,:grafe_5,
                        :grafe_6,:grafe_7,:grafe_8,:grafe_9,:grafe_10,
                        :grafe_11,:grafe_12,:grafe_13,:grafe_14,:grafe_15,
                        :grafe_16,:grafe_17,:grafe_18,:grafe_19,
                        # Year Fixed Effects
                        :year_2016,:year_2017,:year_2018],
    spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
    spec_prodInt=[:padj,:iplan,:inet],
    spec_fixInt=Vector{Symbol}(undef,0),
    spec_fixEff=[:metal_gra,:iss_net_gra],
    spec_wgt= spec_wgt,
    method="ga",ga_itr = 50,ll_start=true)
#
# #
# println("###############################")
# println("Specification 4")
# ### Full Specification - Active Enrollees Only
# filename = "Spec4_$rundate"
# mx_out_1 = MainSpec(df_LA_active,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)

# println("###############################")
# println("Specification 3- Demographic Interaction Terms")
# ### Full Specification
# filename = "Spec3_demint_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,
#             :issfe_1_age_1, :issfe_1_age_2, :issfe_1_fam, :issfe_1_hassub,
#             :issfe_2_age_1, :issfe_2_age_2, :issfe_2_fam, :issfe_2_hassub,
#             :issfe_5_age_1, :issfe_5_age_2, :issfe_5_fam, :issfe_5_hassub,
#             :issfe_6_age_1, :issfe_6_age_2, :issfe_6_fam, :issfe_6_hassub],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,
#                         #Continuous Premium
#                         # :autodp,
#                         # Categorical Premium
#                         :autodp_disc_1,:autodp_disc_2,:autodp_disc_4,:autodp_disc_5,:autodp_disc_6,
#                         #Metal Fixed Effects
#                         :def_mtl_cat,:def_mtl_gld, # Leave Out Bronze
#                         :def_mtl_hdp,:def_mtl_plt,
#                         :def_mtl_slv,:def_mtl_s73,
#                         :def_mtl_s87,:def_mtl_s94,
#                         # Network Fixed Effects
#                         # :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
#                         # :def_issfe_6, :def_issfe_7, # Leave Out LA Care
#                         # :def_netfe_4, :def_netfe_6, # Drop net02, net03
#                         # :def_netfe_9, :def_netfe_12, # Drop net10, net08
#                         #Rating Area Fixed Effect
#                         :grafe_2,:grafe_3,:grafe_4,:grafe_5,
#                         :grafe_6,:grafe_7,:grafe_8,:grafe_9,:grafe_10,
#                         :grafe_11,:grafe_12,:grafe_13,:grafe_14,:grafe_15,
#                         :grafe_16,:grafe_17,:grafe_18,:grafe_19,
#                         # Year Fixed Effects
#                         :year_2016,:year_2017,:year_2018],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)

# println("###############################")
# println("Specification 4, with full demographic interactions")
# ### Full Specification - Active Enrollees Only
# filename = "Spec4_demint_$rundate"
# mx_out_1 = MainSpec(df_LA_active,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,
#             :issfe_1_age_1, :issfe_1_age_2, :issfe_1_fam, :issfe_1_hassub,
#             :issfe_2_age_1, :issfe_2_age_2, :issfe_2_fam, :issfe_2_hassub,
#             :issfe_5_age_1, :issfe_5_age_2, :issfe_5_fam, :issfe_5_hassub,
#             :issfe_6_age_1, :issfe_6_age_2, :issfe_6_fam, :issfe_6_hassub],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)
#
# println("###############################")
# println("Specification 3, with auto in plan choice Robust Check")
# filename = "Spec3_auto_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:newautoplan],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,
#                         #Continuous Premium
#                         # :autodp,
#                         # Categorical Premium
#                         :autodp_disc_1,:autodp_disc_2,:autodp_disc_4,:autodp_disc_5,:autodp_disc_6,
#                         #Metal Fixed Effects
#                         :def_mtl_cat,:def_mtl_gld, # Leave Out Bronze
#                         :def_mtl_hdp,:def_mtl_plt,
#                         :def_mtl_slv,:def_mtl_s73,
#                         :def_mtl_s87,:def_mtl_s94,
#                         # Network Fixed Effects
#                         # :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
#                         # :def_issfe_6, :def_issfe_7, # Leave Out LA Care
#                         # :def_netfe_4, :def_netfe_6, # Drop net02, net03
#                         # :def_netfe_9, :def_netfe_12, # Drop net10, net08
#                         #Rating Area Fixed Effect
#                         :grafe_2,:grafe_3,:grafe_4,:grafe_5,
#                         :grafe_6,:grafe_7,:grafe_8,:grafe_9,:grafe_10,
#                         :grafe_11,:grafe_12,:grafe_13,:grafe_14,:grafe_15,
#                         :grafe_16,:grafe_17,:grafe_18,:grafe_19,
#                         # Year Fixed Effects
#                         :year_2016,:year_2017,:year_2018],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet,:newautoplan],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)
#
#
# println("###############################")
# println("Specification 4, with auto in plan choice Robust Check")
# ### Full Specification
# filename = "Spec4_auto_$rundate"
# mx_out_1 = MainSpec(df_LA_active,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:newautoplan],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet,:newautoplan],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)
#
# println("###############################")
# println("Specification 3, with Kaiser Robust Check")
# filename = "Spec3_kais_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:inet_Kaiser],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,
#                         #Continuous Premium
#                         # :autodp,
#                         # Categorical Premium
#                         :autodp_disc_1,:autodp_disc_2,:autodp_disc_4,:autodp_disc_5,:autodp_disc_6,
#                         #Metal Fixed Effects
#                         :def_mtl_cat,:def_mtl_gld, # Leave Out Bronze
#                         :def_mtl_hdp,:def_mtl_plt,
#                         :def_mtl_slv,:def_mtl_s73,
#                         :def_mtl_s87,:def_mtl_s94,
#                         # Network Fixed Effects
#                         # :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
#                         # :def_issfe_6, :def_issfe_7, # Leave Out LA Care
#                         # :def_netfe_4, :def_netfe_6, # Drop net02, net03
#                         # :def_netfe_9, :def_netfe_12, # Drop net10, net08
#                         #Rating Area Fixed Effect
#                         :grafe_2,:grafe_3,:grafe_4,:grafe_5,
#                         :grafe_6,:grafe_7,:grafe_8,:grafe_9,:grafe_10,
#                         :grafe_11,:grafe_12,:grafe_13,:grafe_14,:grafe_15,
#                         :grafe_16,:grafe_17,:grafe_18,:grafe_19,
#                         # Year Fixed Effects
#                         :year_2016,:year_2017,:year_2018],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)
#
#
# println("###############################")
# println("Specification 4, with Kaiser Robust Check")
# ### Full Specification
# filename = "Spec4_kais_$rundate"
# mx_out_1 = MainSpec(df_LA_active,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet,:inet_Kaiser],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= Vector{Symbol}(undef,0),
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)
#
# println("###############################")
# println("Specification 3, with large threshold cross in attention equation")
# ### Full Specification
# filename = "Spec3_thresh95_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,
#                         #Continuous Premium
#                         # :autodp,
#                         # Categorical Premium
#                         :autodp_disc_1,:autodp_disc_2,:autodp_disc_4,:autodp_disc_5,:autodp_disc_6,
#                         # Large change in the choice set
#                         :choiceset_exceed_thresh_95,:choiceset_below_thresh_05,
#                         #Metal Fixed Effects
#                         :def_mtl_cat,:def_mtl_gld, # Leave Out Bronze
#                         :def_mtl_hdp,:def_mtl_plt,
#                         :def_mtl_slv,:def_mtl_s73,
#                         :def_mtl_s87,:def_mtl_s94,
#                         # Network Fixed Effects
#                         # :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
#                         # :def_issfe_6, :def_issfe_7, # Leave Out LA Care
#                         # :def_netfe_4, :def_netfe_6, # Drop net02, net03
#                         # :def_netfe_9, :def_netfe_12, # Drop net10, net08
#                         #Rating Area Fixed Effect
#                         :grafe_2,:grafe_3,:grafe_4,:grafe_5,
#                         :grafe_6,:grafe_7,:grafe_8,:grafe_9,:grafe_10,
#                         :grafe_11,:grafe_12,:grafe_13,:grafe_14,:grafe_15,
#                         :grafe_16,:grafe_17,:grafe_18,:grafe_19,
#                         # Year Fixed Effects
#                         :year_2016,:year_2017,:year_2018],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)
#
# println("###############################")
# println("Specification 3, with large threshold cross in attention equation")
# ### Full Specification
# filename = "Spec3_thresh98_$rundate"
# mx_out_1 = MainSpec(df_LA,filename,
#     haltonDim = mixed_draws,
#     spec_per = spec_per,
#     spec_prd = spec_prd,
#     spec_ch = spec_ch,
#     spec_ch_last = spec_ch_last,
#     spec_prodchr = [:padj,:iplan,:inet],
#     spec_prodchr_0= [:issfe_1, :issfe_2, :issfe_5, :issfe_6],
#     spec_inertchr= [:constant,:agefe_1,:agefe_2,:fam,:hassub,
#                         #Continuous Premium
#                         # :autodp,
#                         # Categorical Premium
#                         :autodp_disc_1,:autodp_disc_2,:autodp_disc_4,:autodp_disc_5,:autodp_disc_6,
#                         # Large change in the choice set
#                         :choiceset_exceed_thresh_98,:choiceset_below_thresh_02,
#                         #Metal Fixed Effects
#                         :def_mtl_cat,:def_mtl_gld, # Leave Out Bronze
#                         :def_mtl_hdp,:def_mtl_plt,
#                         :def_mtl_slv,:def_mtl_s73,
#                         :def_mtl_s87,:def_mtl_s94,
#                         # Network Fixed Effects
#                         # :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
#                         # :def_issfe_6, :def_issfe_7, # Leave Out LA Care
#                         # :def_netfe_4, :def_netfe_6, # Drop net02, net03
#                         # :def_netfe_9, :def_netfe_12, # Drop net10, net08
#                         #Rating Area Fixed Effect
#                         :grafe_2,:grafe_3,:grafe_4,:grafe_5,
#                         :grafe_6,:grafe_7,:grafe_8,:grafe_9,:grafe_10,
#                         :grafe_11,:grafe_12,:grafe_13,:grafe_14,:grafe_15,
#                         :grafe_16,:grafe_17,:grafe_18,:grafe_19,
#                         # Year Fixed Effects
#                         :year_2016,:year_2017,:year_2018],
#     spec_demR=[:agefe_1,:agefe_2,:fam,:hassub],
#     spec_prodInt=[:padj,:iplan,:inet],
#     spec_fixInt=Vector{Symbol}(undef,0),
#     spec_fixEff=[:metal_gra,:iss_net_gra],
#     spec_wgt= spec_wgt,
#     method="ga",ga_itr = 50,ll_start=true)

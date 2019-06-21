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
### Add Price/Default Interaction
df_LA[:padj_defpadj] = df_LA[:padj].*df_LA[:def_padj]

println("Data Loaded")


rundate = "2019-05-25"
file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/ML_spec6_$rundate.jld2"
@load file p_est spec_Dict fval

# ## Full Model
# Structure the data
c = ChoiceData(df_LA;
    per = [:hh_year_id],
    prd = [:product],
    ch = [:choice],
    ch_last = [:iplan],
    prodchr = [:padj,:iplan,:inet,:iiss,
    :issfe_1, :issfe_2, :issfe_5, :issfe_6,
    :issfe_8, :issfe_9, # Leave Out LA Care
    :netfe_2, :netfe_3, :netfe_4, :netfe_7,
    :netfe_11, :netfe_12, :netfe_13, :netfe_15],
    prodchr_0=[:issfe_1, :issfe_2, :issfe_5, :issfe_6],
    inertchr=[:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,:def_padj,
                    :def_mtl_brz,:def_mtl_cat,:def_mtl_gld,
                    :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
                    :def_mtl_s87,:def_mtl_s94],
    # inertchr=Vector{Symbol}(undef,0),
    demR =[:agefe_1,:agefe_2,:fam,:hassub],
    prodInt=[:padj,:iplan,:inet,:iiss],
    fixEff=[:metal],
    wgt=[:constant])

# Fit into model
m = InsuranceLogit(c,500)

ReturnPercBase, ReturnPercObs = predict_switching(m,p_est)
println(ReturnPercObs)
println(ReturnPercBase)
ReturnFA, ReturnPercObs = predict_switching(m,p_est,fullAtt=true)
println(ReturnFA)
ReturnNoCont, ReturnPercObs = predict_switching(m,p_est,noCont=true)
println(ReturnNoCont)
ReturnNoHass, ReturnPercObs = predict_switching(m,p_est,noHass=true)
println(ReturnNoHass)
ReturnNone, ReturnPercObs = predict_switching(m,p_est,fullAtt=true,noHass=true)
println(ReturnNone)
ReturnNone, ReturnPercObs = predict_switching(m,p_est,fullAtt=true,noHass=true,noCont=true)
println(ReturnNone)



##### Check Active Relationship ####

active_obs, active_pred = activePredict(m,p_est,df_LA)

out1 = DataFrame(active_obs=active_obs,active_pred=active_pred)
file1 = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/active_$rundate.csv"
CSV.write(file1,out1)

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
df_LA = df
 # df_LA = df[df[:gra].==16,:]
# df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")
println("Data Loaded")


rundate = "2019-08-02"
spec = "Spec1b_"
file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/$spec$rundate.jld2"
@load file p_est spec_Dict fval

# ## Full Model
# Structure the data
c = ChoiceData(df_LA;
    per = [:hh_id],
    # per = spec_Dict["per"],
    prd = spec_Dict["prd"],
    ch = spec_Dict["ch"],
    ch_last = spec_Dict["ch_last"],
    prodchr = spec_Dict["prodchr"],
    prodchr_0=spec_Dict["prodchr_0"],
    inertchr=spec_Dict["inertchr"],
    # inertchr=Vector{Symbol}(undef,0),
    demR =spec_Dict["demR"],
    prodInt=spec_Dict["prodInt"],
    fixEff=spec_Dict["fixEff"],
    wgt=spec_Dict["wgt"])

# Fit into model
m = InsuranceLogit(c,spec_Dict["haltonDim"])
if m.parLength[:All]!=length(p_est)
    error("WARNING: Specification Error!")
end

ll = log_likelihood(m,p_est)
println(ll)

parBase = parDict(m,p_est)
individual_values!(m,parBase)
individual_shares(m,parBase)
# grad = Vector{Float64}(undef,length(p_est))


ReturnPercBase, ReturnPercObs = predict_switching(m,p_est,spec_Dict)
println(ReturnPercObs)
println(ReturnPercBase)
ReturnFA, ReturnPercObs = predict_switching(m,p_est,spec_Dict,fullAtt=true)
println(ReturnFA)
ReturnNoHass, ReturnPercObs = predict_switching(m,p_est,spec_Dict,noHass=true)
println(ReturnNoHass)
ReturnNoCont, ReturnPercObs = predict_switching(m,p_est,spec_Dict,noCont=true)
println(ReturnNoCont)
ReturnNone, ReturnPercObs = predict_switching(m,p_est,spec_Dict,fullAtt=true,noHass=true)
println(ReturnNone)
ReturnNone, ReturnPercObs = predict_switching(m,p_est,spec_Dict,fullAtt=true,noHass=true,noCont=true)
println(ReturnNone)

#### Average Willingness to Pay ####


βMat = coeff_values(m,parBase)
wtp_iplan = -100*(βMat[:,2]./βMat[:,1])
wtp_inet = -100*(βMat[:,3]./βMat[:,1])
wtp_iiss = -100*(βMat[:,4]./βMat[:,1])
wtp_cont = -100*((βMat[:,2]+βMat[:,3] .+ βMat[:,4])./βMat[:,1])

println("Plan Level")
println(mean(wtp_iplan))
println(std(wtp_iplan))
println("Network Level")
println(mean(wtp_inet))
println(std(wtp_inet))
println("Issuer Level")
println(mean(wtp_iiss))
println(std(wtp_iiss))
println("Total Continuity")
println(mean(wtp_cont))
println(std(wtp_cont))


##### Check Active Relationship ####

switchers,returning, active_obs, active_pred = activePredict(m,p_est,df_LA)

out1 = DataFrame(switch=switchers,returning=returning,active_obs=active_obs,active_pred=active_pred)
file1 = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/active_$rundate.csv"
CSV.write(file1,out1)

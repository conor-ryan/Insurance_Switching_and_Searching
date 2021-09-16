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
df_LA = df
# df_LA = df_active
 # df_LA = df[df[:gra].==8,:]
 # df_LA = df[df[:gra].==16,:]
# df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")
println("Data Loaded")


rundate = "2021-08-27"
spec = "Spec3_"
file = "$home/Output/Estimation_Results/$spec$rundate.jld2"
@load file p_est spec_Dict fval

println("Running Spec: $spec$rundate")

useActiveVar=false
if spec=="Spec4_"
    useActiveVar=true
end

# ## Full Model
# Structure the data
c = ChoiceData(df_LA;
    # per = [:hh_id],
    per = spec_Dict["per"],
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

# numPar = length(p_est)
# par = parDict(m,p_est)
# Pop = length(par.ω_i)
# ll = fval*Pop
# BIC = log(Pop)*(numPar+1) - 2*ll
# println("Number of Parameters: $numPar, Log-Likelihood: $ll, BIC: $BIC")
#
#
#
# ll = log_likelihood(m,p_est)
# println(ll*Pop)
#
# grad = Vector{Float64}(undef,length(p_est))

# println("Base")
# resBase = predict_switching(m,p_est,spec_Dict,useActiveVar=useActiveVar)
# println("Full Attention")
# resFA = predict_switching(m,p_est,spec_Dict,fullAtt=true,useActiveVar=useActiveVar)
# println("No Hassle Costs")
# resNoHass = predict_switching(m,p_est,spec_Dict,noHass=true,useActiveVar=useActiveVar)
# println("No Taste for Continuity")
# resNoCont = predict_switching(m,p_est,spec_Dict,noCont=true,useActiveVar=useActiveVar)
# println("Only Inattention")
# resOnlyAtt = predict_switching(m,p_est,spec_Dict,noHass=true,noCont=true,useActiveVar=useActiveVar)
# println("Only Hassle Costs")
# resOnlyHass = predict_switching(m,p_est,spec_Dict,fullAtt=true,noCont=true,useActiveVar=useActiveVar)
# println("Only Tastes for Continuity")
# resOnlyCont = predict_switching(m,p_est,spec_Dict,fullAtt=true,noHass=true,useActiveVar=useActiveVar)
# println("No Switching Costs")
# resNone = predict_switching(m,p_est,spec_Dict,fullAtt=true,noHass=true,noCont=true,useActiveVar=useActiveVar)
#
# #### Average Willingness to Pay ####
parBase = parDict(m,p_est)
individual_values!(m,parBase)
individual_shares(m,parBase)

#### Save predicted Choices #####
s_ij = parBase.s_hat

out = DataFrame(Person = df_LA[!,:hh_id], Product = df_LA[!,:product], Year = df_LA[!,:year],
                    Metal = df_LA[!,:metal], Issuer = df_LA[!,:issuername],
                    s_pred = parBase.s_hat)
file1 = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/predictedChoices_$spec$rundate.csv"
CSV.write(file1,out)



#
#
# βMat = coeff_values(m,parBase)
# wtp_iplan = -100*(βMat[:,2]./βMat[:,1])
# wtp_inet = -100*(βMat[:,3]./βMat[:,1])
# wtp_iiss = -100*(βMat[:,4]./βMat[:,1])
# wtp_cont = -100*((βMat[:,2]+βMat[:,3])./βMat[:,1])


# alpha_long = parBase.β_0[1] .+ parBase.β[1,:]'*demoRaw(m.data)
# price_long = prodchars(m.data)[1,:]
# elas = Vector{Float64}(undef,length(alpha_long))
# for i in 1:length(elas)
#     elas[i] = alpha_long[i]*price_long[i]*(1 - parBase.s_hat[i])
# end

# println("Plan Level")
# println(mean(wtp_iplan))
# println(std(wtp_iplan))
# println("Network Level")
# println(mean(wtp_inet))
# println(std(wtp_inet))
# println("Issuer Level")
# println(mean(wtp_iiss))
# println(std(wtp_iiss))
# println("Total Continuity")
# println(mean(wtp_cont))
# println(std(wtp_cont))

# ## Marginal Effects
# println("### COMPUTE MARGINAL EFFECTS ### ")
# println(round.(100 .*marginalEffects(m,p_est),digits=2))
#
# dME = 100 .*meDeriv(m,p_est)
# Var, se1, se2,t_stat, stars = res_process(m,p_est)
# Σ_coeff = Diagonal(se1[1:m.parLength[:I]])
# Σ_me = dME*Var[1:m.parLength[:I],1:m.parLength[:I]]*dME'
# σ_me = sqrt.(diag(Σ_me))
# println("### MARGINAL EFFECT STANDARD ERRORS ###")
# println(round.(σ_me,digits=2))

# ##### Demographic Buckets #####
# β0 = parBase.β_0[1:3]
# β = parBase.β[1:3,:]
# Z =[0 0 0 1;
#     1 0 0 1;
#     0 1 0 1;
#     0 0 1 1;
#     1 0 1 1;
#     0 1 1 1;
#     0 0 0 0;
#     1 0 0 0;
#     0 1 0 0;
#     0 0 1 0;
#     1 0 1 0;
#     0 1 1 0;]
#
# βz = β0 .+ β*Z'
# J,K = size(βz)
# WTP = zeros(K,J-1)
# for j in 1:(J-1), k in 1:K
#     WTP[k,j] = -100*βz[1+j,k]/βz[1,k]
# end
# out = DataFrame(WTP)
#
# file1 = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/wtp_$spec$rundate.csv"
# CSV.write(file1,out)
#
# ##### Check Active Relationship ####
# par = parDict(m,p_est)
# individual_values!(m,par)
# individual_shares(m,par)
# activePredict(m,par,df_LA,spec,rundate)

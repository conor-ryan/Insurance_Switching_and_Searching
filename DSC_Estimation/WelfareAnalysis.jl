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

#  df_LA = df[df[:gra].==16,:]
# df_active = 0.0
# df = 0.0

# df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")
println("Data Loaded")


rundate = "2020-12-04"
spec = "Spec3_"
file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/$spec$rundate.jld2"
@load file p_est spec_Dict fval



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

Istart = rand(m.parLength[:I])/10 .-.05
βstart = rand(m.parLength[:β])/10 .-.05
γstart = rand(m.parLength[:γ]*length(m.data._prodInteract))/10 .-.05
σstart = rand(m.parLength[:σ])/10 .- .05
FEstart = rand(m.parLength[:FE])/100 .-.005

p_est = vcat(Istart,βstart,γstart,σstart,FEstart)

numPar = length(p_est)
par = parDict(m,p_est)
Pop = length(par.ω_i)
ll = fval*Pop
BIC = log(Pop)*(numPar+1) - 2*ll
println("Number of Parameters: $numPar, Log-Likelihood: $ll, BIC: $BIC")



ll = log_likelihood(m,p_est)
println(ll*Pop)

grad = Vector{Float64}(undef,length(p_est))


function all_shares(d::InsuranceLogit,p::parDict{T}) where T
    # Store Parameters
    μ_ij_large = p.μ_ij
    s_ij = similar(μ_ij_large)
    ω_large = p.ω_i
    iplan_large = choice_last(d.data)[:]
    y_large = choice(d.data)[:]

    for (ind,idxitr) in d.data._personDict
        # println(ind)
        u = μ_ij_large[:,idxitr]
        ω = ω_large[d.data._searchDict[ind]]
        iplan = iplan_large[idxitr]
        dict = d.data._personYearDict[ind]
        yr_next = findYearInd(d.data._personYearDict[ind])
        y = y_large[idxitr]
        y = findall(y.==1)
        s_hat,s_cond,ll,expsum =calc_shares_mat(u,ω,iplan,y,yr_next)
        s_ij[:,idxitr] = transpose(s_hat)
    end
    return s_ij
end

function remove_switching_pars(m::InsuranceLogit,p_vec::Vector{Float64},spec::Dict{String,Any};
    fullAtt::Bool=false,
    noCont::Bool=false,
    noHass::Bool=false)

    p_est = copy(p_vec)

    Ilength = m.parLength[:I]

    parBase = parDict(m,p_est)

    if noCont
        β_ind = inlist(spec["prodchr"],[:inet,:iiss])
        parBase.β[β_ind,:].=0.0
        parBase.β_0[β_ind].= 0.0
    end
    if noHass
        β_ind = findall(spec["prodchr"].==:iplan)
        parBase.β[β_ind,:].=0.0
        parBase.β_0[β_ind].= 0.0
    end

    individual_values!(m,parBase)

    if fullAtt
        parBase.ω_i[:] .= 1.0
    end

    individual_shares(m,parBase)

    return parBase
end

#### Consumer Alpha ####
parBase = parDict(m,p_est)
β0 = parBase.β_0[1:3]
β = parBase.β[1:3,:]
Z = demoRaw(c)
α = (β0 .+ β*Z)[1,:]
println("Median alpha is $(median(α))")
#### Switching Cost Neutral WTP ####
ret_index = returning_index(m,parBase)
parBase = nothing

parNeutral = remove_switching_pars(m,p_est,spec_Dict,noHass=true)
μ_ij = transpose(parNeutral.μ_ij)

WTP_insurance = similar(μ_ij)

for i in 1:size(μ_ij,1), j in 1:size(μ_ij,2)
    WTP_insurance[i,j] = -μ_ij[i,j]/α[j]
end
parNeutral = nothing

WTP_insurance = WTP_insurance[ret_index,:]


#### Test Highest and Lowest WTP ####
wtp_max = maximum(WTP_insurance,dims=2)
wtp_min = minimum(WTP_insurance,dims=2)
at_stake = wtp_max - wtp_min
println("Average at stake WTP: $(mean(at_stake))")

#### Choice Probabilities #####
function return_choices(m::InsuranceLogit,p_vec::Vector{Float64},spec::Dict{String,Any},ret_index::Vector{Int},WTP::Matrix{Float64};
    fullAtt::Bool=false,
    noCont::Bool=false,
    noHass::Bool=false)
    println("Set Parameters")
    par = remove_switching_pars(m,p_est,spec_Dict,fullAtt=fullAtt,noCont=noCont,noHass=noHass)
    # individual_values!(m,par)
    # individual_shares(m,par)
    shares = all_shares(m,par)
    shares = transpose(shares)
    shares = shares[ret_index,:]
    println("Compute Mean")
    avg_wtp = mean(shares.*WTP)
    return avg_wtp
end



wtp_base = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance)
println("Base WTP: $wtp_base")
wtp_noHass = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,noHass=true)
println("WTP, no Hass: $wtp_noHass")
wtp_noCont = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,noCont=true)
println("WTP, no Continuity: $wtp_noCont")
wtp_fullAtten = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,fullAtt=true)
println("WTP, full attention: $wtp_fullAtten")
wtp_fullAtten_noHass = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,noHass=true,fullAtt=true)
println("WTP, full attention & no hassle costs: $wtp_fullAtten_noHass")
wtp_fullAtten_noCont = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,noCont=true,fullAtt=true)
println("WTP, full attention & no continuity: $wtp_fullAtten_noCont")
wtp_noHass_noCont = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,noCont=true,noHass=true)
println("WTP, no hassle & continuity: $wtp_noHass_noCont")

using JLD2
using CSV
using Random
using Dates
using LinearAlgebra
using Statistics
using Distributions

Random.seed!(123)

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

function draw_WTP(d::InsuranceLogit,p::parDict{T},
    eps_draws::Matrix{Float64},search_draws::Vector{Float64},WTP::Matrix{Float64}) where T
    # Store Parameters
    μ_ij_large = log.(p.μ_ij)
    ω_large = p.ω_i
    WTP_byperson = Vector{Float64}(undef,length(p.ω_i))
    # println("Mean Atten: $(mean(ω_large))")
    iplan_large = choice_last(d.data)[:]
    y_large = choice(d.data)[:]

    for (ind,idxitr) in d.data._personDict
        # println(ind)
        u = μ_ij_large[:,idxitr] + eps_draws[:,idxitr]
        WTP_per = WTP[:,idxitr]
        ω_index = d.data._searchDict[ind]
        ω = ω_large[ω_index]
        ω_draws = search_draws[ω_index]
        iplan = iplan_large[idxitr]
        dict = d.data._personYearDict[ind]
        y = y_large[idxitr]
        y = findall(y.==1)
        years = sort(Int.(keys(d.data._personYearDict[ind])))
        yr_ind = 0
        for year in years
            idx_year = d.data._personYearDict[ind][year]
            yr_ind+=1
            ω_year = ω[yr_ind]
            util = u[:,idx_year]
            max_util = maximum(util,dims=2)
            share_cond = zeros(size(util))


            ω_year_draw = ω_draws[yr_ind]
            if ω_year_draw<ω_year
                for n in 1:size(share_cond,1)
                    share_cond[n,:] = Float64.(util[n,:].==max_util[n])
                end
            else
                for n in 1:size(share_cond,1)
                    share_cond[n,:] = iplan[idx_year]
                end
            end
            mean_WTP = mean(sum(WTP_per[:,idx_year].*share_cond,dims=2))
            WTP_byperson[ω_index[yr_ind]] = mean_WTP
        end
    end
    return WTP_byperson
end

function remove_switching_pars(m::InsuranceLogit,p_vec::Vector{Float64},spec::Dict{String,Any};
    fullAtt::Bool=false,
    noCont::Bool=false,
    noHass::Bool=false,
    useActiveVar::Bool=false)

    p_est = copy(p_vec)

    Ilength = m.parLength[:I]

    par = parDict(m,p_est)
    if useActiveVar
        par.ω_i[:] = m.data.active[:]
    end
    # test = mean(par.ω_i)
    # println("Test Inattention 1: $test")

    if noCont
        β_ind = inlist(spec["prodchr"],[:inet,:iiss])
        par.β[β_ind,:].=0.0
        par.β_0[β_ind].= 0.0
    end
    if noHass
        β_ind = findall(spec["prodchr"].==:iplan)
        par.β[β_ind,:].=0.0
        par.β_0[β_ind].= 0.0
    end

    individual_values!(m,par)

    if fullAtt
        par.ω_i[:] .= 1.0
    end

    # test = mean(par.ω_i)
    # println("Test Inattention 2: $test")

    individual_shares(m,par)

    return par
end

function return_choices(m::InsuranceLogit,p_vec::Vector{Float64},spec::Dict{String,Any},
    ret_index::Vector{Int},WTP::Matrix{Float64},eps_draws::Matrix{Float64},search_draws::Vector{Float64};
    fullAtt::Bool=false,
    noCont::Bool=false,
    noHass::Bool=false,
    useActiveVar::Bool=false)
    println("Set Parameters")
    par = remove_switching_pars(m,p_est,spec_Dict,fullAtt=fullAtt,noCont=noCont,noHass=noHass,
    useActiveVar=useActiveVar)
    # individual_values!(m,par)
    # individual_shares(m,par)
    wtp = draw_WTP(m,par,eps_draws,search_draws,WTP)
    avg_wtp = mean(wtp[ret_index])
    return avg_wtp
end


# Load the Data
include("load.jl")
df_LA = df

# df_LA = df[df[:gra].==10,:]
df_active = 0.0
df = 0.0

# df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")
println("Data Loaded")

# rundate = "2020-12-14"
rundate = "2021-01-31"
spec = "Spec3_"
file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/$spec$rundate.jld2"
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

numPar = length(p_est)
par = parDict(m,p_est)
Pop = length(par.ω_i)
ll = fval*Pop
BIC = log(Pop)*(numPar+1) - 2*ll
println("Number of Parameters: $numPar, Log-Likelihood: $ll, BIC: $BIC")



ll = log_likelihood(m,p_est)
println(ll*Pop)

grad = Vector{Float64}(undef,length(p_est))

#### Consumer Alpha ####
parBase = parDict(m,p_est)
β0 = parBase.β_0[1:3]
β = parBase.β[1:3,:]
Z = demoRaw(c)
α = (β0 .+ β*Z)[1,:]
println("Median alpha is $(median(α))")

# individual_values!(m,parBase)
# μ_ij = log.(parBase.μ_ij)
# println("Mean of all utility values: $(mean(μ_ij))")
# println("Median of all utility values: $(median(μ_ij))")
# println("Maximum of all utility values: $(maximum(μ_ij))")
# println("Minimum of all utility values: $(minimum(μ_ij))")
#
# individual_shares(m,parBase)
# weighted_utility = parBase.s_hat.*μ_ij
# println("Mean of weighted utility?: $(mean(weighted_utility))")

#### Epsilon Draws ####
t1ev = Gumbel(0,1)
eps_draws = rand(t1ev,size(parBase.μ_ij))
search_draws = rand(length(parBase.ω_i))

#### Switching Cost Neutral WTP ####
ret_index = returning_peryear_index(m,parBase)
parBase = nothing

parNeutral = remove_switching_pars(m,p_est,spec_Dict,noHass=true)


μ_ij = log.(parNeutral.μ_ij)
println("Mean of all utility values: $(mean(μ_ij))")

WTP_insurance = similar(μ_ij)

for i in 1:size(μ_ij,1), j in 1:size(μ_ij,2)
    WTP_insurance[i,j] = -μ_ij[i,j]/α[j]
end
parNeutral = nothing


#### Test Highest and Lowest WTP ####
wtp_max = maximum(WTP_insurance,dims=1)
wtp_min = minimum(WTP_insurance,dims=1)
at_stake = wtp_max - wtp_min
println("Average at stake WTP: $(mean(at_stake))")

### Choice Probabilities #####
wtp_base = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,useActiveVar=useActiveVar)
println("Base WTP: $wtp_base")
wtp_noHass = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,noHass=true,useActiveVar=useActiveVar)
println("WTP, no Hass: $wtp_noHass")
wtp_noCont = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,noCont=true,useActiveVar=useActiveVar)
println("WTP, no Continuity: $wtp_noCont")
wtp_fullAtten = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,fullAtt=true,useActiveVar=useActiveVar)
println("WTP, full attention: $wtp_fullAtten")
wtp_fullAtten_noHass = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,noHass=true,fullAtt=true,useActiveVar=useActiveVar)
println("WTP, full attention & no hassle costs: $wtp_fullAtten_noHass")
wtp_fullAtten_noCont = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,noCont=true,fullAtt=true,useActiveVar=useActiveVar)
println("WTP, full attention & no continuity: $wtp_fullAtten_noCont")
wtp_noHass_noCont = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,noCont=true,noHass=true,useActiveVar=useActiveVar)
println("WTP, no hassle & continuity: $wtp_noHass_noCont")
wtp_noCosts = return_choices(m,p_est,spec_Dict,ret_index,WTP_insurance,eps_draws,search_draws,noCont=true,noHass=true,fullAtt=true,useActiveVar=useActiveVar)
println("WTP, no switching costs: $wtp_noCosts")

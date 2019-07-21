using JLD2
using CSV
using Random
using Dates
using LinearAlgebra
using Statistics
using BenchmarkTools
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
# 1% Sample of LA
df_LA = df
# df = df[df[:samp1].==1,:]
# df_LA = df[df[:gra].==16,:]


### Add Anthem Fixed Effect
# df_LA[:issfe_1] = Int.(df_LA[:issuername].=="Anthem")
### Add Price/Default Interaction
# df_LA[:padj_defpadj] = df_LA[:padj].*df_LA[:def_padj]

println("Data Loaded")
# df = 0.0

# Structure the data
c = ChoiceData(df_LA;
    per = [:hh_id],
    prd = [:product],
    ch = [:choice],
    ch_last = [:iplan],
    prodchr =  [:padj,:iplan,:inet,:iiss,
    :issfe_1, :issfe_2, :issfe_3, :issfe_4,
    :issfe_6, :issfe_7, # Leave Out LA Care
    :netfe_2, :netfe_3, :netfe_4, :netfe_6,
    :netfe_8, :netfe_9, :netfe_10, :netfe_12],
    # prodchr_0=Vector{Symbol}(undef,0),
    prodchr_0=[:issfe_1, :issfe_2, :issfe_3, :issfe_4],
    inertchr=[:constant,:agefe_1,:agefe_2,:fam,:hassub,:dprem,
                        #Metal Fixed Effects
                        :def_mtl_brz,:def_mtl_cat,:def_mtl_gld, # Leave Out Silver
                        :def_mtl_hdp,:def_mtl_plt,:def_mtl_s73,
                        :def_mtl_s87,:def_mtl_s94,
                        # Network Fixed Effects
                        :def_issfe_1, :def_issfe_2, :def_issfe_3, :def_issfe_4,
                        :def_issfe_6, :def_issfe_7, # Leave Out LA Care
                        :def_netfe_2, :def_netfe_4, :def_netfe_6, # Drop net03
                        :def_netfe_8, :def_netfe_9, :def_netfe_12, # Drop net10
                        # Year Fixed Effects
                        :year_2015,:year_2016,:year_2017,:year_2018],
    demR =[:agefe_1,:agefe_2,:fam,:hassub],
    prodInt=[:padj,:iplan,:inet,:iiss],
    fixEff=[:metal],
    wgt=[:constant])

# Fit into model
m = InsuranceLogit(c,50)
println("Data Loaded")

#γ0start = rand(1)-.5


Istart = rand(m.parLength[:I])/10 .-.05
βstart = rand(m.parLength[:β])/10 .-.05
γstart = rand(m.parLength[:γ]*length(m.data._prodInteract))/10 .-.05
σstart = rand(m.parLength[:σ])/10 .- .05
FEstart = rand(m.parLength[:FE])/100 .-.005

p0 = vcat(Istart,βstart,γstart,σstart,FEstart)

par = parDict(m,p0)

individual_values!(m,par)
individual_shares(m,par)



println("Compute Gradient")
grad = Vector{Float64}(undef,length(p0))
hess = Matrix{Float64}(undef,length(p0),length(p0))
ll = log_likelihood!(hess,grad,m,p0)

# app = iterate(eachperson(m.data),7)[1]
# ll = test_grad!(hess,grad,app,m,p0)


# # #
# # #
f_obj(x) = log_likelihood(m,x)
# f_obj(x) = test_grad(app,m,x)
grad_1 = Vector{Float64}(undef,length(p0))
hess_1 = Matrix{Float64}(undef,length(p0),length(p0))
fval_old = f_obj(p0)
println(fval_old-ll)
# # # # # #
println("Grad")
ForwardDiff.gradient!(grad_1,f_obj, p0)#, cfg)
println(maximum(abs.(grad_1-grad)))
grad_ind = findall(abs.(grad_1-grad).>1e-12)
println("Hessian")
# cfg = ForwardDiff.HessianConfig(f_obj, p0, ForwardDiff.Chunk{3}())
ForwardDiff.hessian!(hess_1,f_obj, p0)
println(maximum(abs.(hess_1-hess)))



using Profile
Profile.init(n=10^8,delay=.001)
Profile.clear()
#Juno.@profile add_obs_mat!(hess,grad,hess_obs,grad_obs,Pop)
# Juno.@profile log_likelihood!(thD_2,hess_2,grad_2,m,par0)
Juno.@profile log_likelihood!(hess,grad,m,p0)
Juno.profiletree()
Juno.profiler()


#
par0 = parDict(m,p0)
individual_values!(m,par0)
individual_shares(m,par0)
app = iterate(eachperson(m.data),5)[1]
#
#
# ll = log_likelihood(m,par0)

println("Gradient Test")
# p0 = est[1]
par0 = parDict(m,p0)
grad_1 = Vector{Float64}(undef,length(p0))
grad_2 = Vector{Float64}(undef,length(p0))
hess_2 = Matrix{Float64}(undef,length(p0),length(p0))
ll = log_likelihood!(grad_1,m,p0)

m = InsuranceLogit(c,50)
ll = log_likelihood!(hess_2,grad_1,m,p0)
@time log_likelihood!(hess_2,grad_1,m,p0)
@time log_likelihood!(hess_2,grad_1,m,p0)






p_ga, fval = gradient_ascent(m,p0,grad_tol = 1e-5,x_tol=1e-4,max_itr=100)
p_est,fval,flag = newton_raphson_ll(m,p_ga)


est = newton_raphson_ll(m,x_start)
p_est, fval = est

Pop = sum(weight(m.data).*choice(m.data))
println(fval*Pop)



par0 = parDict(m,p_est)
individual_values!(m,par0)
individual_shares(m,par0)


using Profile
Profile.init(n=10^8,delay=.001)
Profile.clear()
#Juno.@profile add_obs_mat!(hess,grad,hess_obs,grad_obs,Pop)
# Juno.@profile log_likelihood!(thD_2,hess_2,grad_2,m,par0)
Juno.@profile res =  log_likelihood!(hess_2,grad_2,m,p0)
Juno.profiletree()
Juno.profiler()



#### Results Implications ####
file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/ML_spec6_2019-05-25.jld2"
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
    demR =[:agefe_1,:agefe_2,:fam,:hassub],
    prodInt=[:padj,:iplan,:inet,:iiss],
    fixEff=[:metal],
    wgt=[:constant])

# Fit into model
m = InsuranceLogit(c,500)

ReturnPercBase, ReturnPercObs = predict_switching(m,p_est)








Ilength = m.parLength[:I]

cont_pars = vcat(m.parLength[:I] .+ (2:4),(Ilength.+ m.parLength[:β]).+ vcat((2:4), 4 .+ (2:4),8 .+ (2:4),12 .+ (2:4)))

println(p_est[1:Ilength])
println(p_est[cont_pars])

#### Baseline ####
parBase = parDict(m,p_est)
individual_values!(m,parBase)
individual_shares(m,parBase)

inertPlan = choice_last(m.data)
obsPlan = choice(m.data)

#### Calculate Totals of Returning Individuals ####
All_Return = [0.0]
All_Stay = [0.0]
All_Stay_Obs = [0.0]
for i in m.data._personIDs
    idx = m.data._personDict[i]
    returning = sum(inertPlan[idx])
    if returning==0.0
        continue
    end
    All_Return[:] = All_Return[:] .+ 1
    retplan = findall(inertPlan[idx].>0)
    All_Stay[:] =All_Stay[:] .+ parBase.s_hat[idx[retplan]]
    obs = findall(obsPlan[idx].>0)
    if retplan==obs
        All_Stay_Obs[:] = All_Stay_Obs[:] .+ 1
    end
end

ReturnPercBase = All_Stay[1]/All_Return[1]
ReturnPercObs = All_Stay_Obs[1]/All_Return[1]


## Full Consideration
parBase.ω_i.=1.0
individual_shares(m,parBase)

All_Return = [0.0]
All_Stay = [0.0]
All_Stay_Obs = [0.0]
for i in m.data._personIDs
    idx = m.data._personDict[i]
    returning = sum(inertPlan[idx])
    if returning==0.0
        continue
    end
    All_Return[:] = All_Return[:] .+ 1
    retplan = findall(inertPlan[idx].>0)
    All_Stay[:] =All_Stay[:] .+ parBase.s_hat[idx[retplan]]
    obs = findall(obsPlan[idx].>0)
    if retplan==obs
        All_Stay_Obs[:] = All_Stay_Obs[:] .+ 1
    end
end

ReturnFuAtt = All_Stay[1]/All_Return[1]



#### No Switching Cost ####
p_noCont = copy(p_est)
p_noCont[cont_pars] .= 0.0
parNoCt = parDict(m,p_noCont)
individual_values!(m,parNoCt)
individual_shares(m,parNoCt)

All_Return = [0.0]
All_Stay = [0.0]
All_Stay_Obs = [0.0]
for i in m.data._personIDs
    idx = m.data._personDict[i]
    returning = sum(inertPlan[idx])
    if returning==0.0
        continue
    end
    All_Return[:] = All_Return[:] .+ 1
    retplan = findall(inertPlan[idx].>0)
    All_Stay[:] =All_Stay[:] .+ parNoCt.s_hat[idx[retplan]]
    obs = findall(obsPlan[idx].>0)
    if retplan==obs
        All_Stay_Obs[:] = All_Stay_Obs[:] .+ 1
    end
end

ReturnNoCont = All_Stay[1]/All_Return[1]



## No Inertia
parNoCt.ω_i.=1.0
individual_shares(m,parNoCt)

All_Return = [0.0]
All_Stay = [0.0]
All_Stay_Obs = [0.0]
for i in m.data._personIDs
    idx = m.data._personDict[i]
    returning = sum(inertPlan[idx])
    if returning==0.0
        continue
    end
    All_Return[:] = All_Return[:] .+ 1
    retplan = findall(inertPlan[idx].>0)
    All_Stay[:] =All_Stay[:] .+ parNoCt.s_hat[idx[retplan]]
    obs = findall(obsPlan[idx].>0)
    if retplan==obs
        All_Stay_Obs[:] = All_Stay_Obs[:] .+ 1
    end
end

ReturnNone = All_Stay[1]/All_Return[1]


println(ReturnPercObs)
println(ReturnPercBase)
println(ReturnNoIn)
println(ReturnNoCont)
println(ReturnNone)

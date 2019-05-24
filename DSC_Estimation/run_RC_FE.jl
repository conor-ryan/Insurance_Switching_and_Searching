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
# df = 0.0

# Structre the data
c = ChoiceData(df_LA;
    per = [:hh_year_id],
    prd = [:product],
    ch = [:choice],
    ch_last = [:iplan],
    prodchr = [:padj,:padj_defpadj],
    prodchr_0= Vector{Symbol}(undef,0),
    inertchr=Vector{Symbol}(undef,0),
    demR =Vector{Symbol}(undef,0),
    prodInt=Vector{Symbol}(undef,0),
    fixEff=[:def_metal,:def_netfe],
    # fixEff=Vector{Symbol}(undef,0),
    fixInt=[:metal,:netname],
    wgt=[:constant])

# Fit into model
m = InsuranceLogit(c,1)
println("Data Loaded")

#γ0start = rand(1)-.5


Istart = rand(m.parLength[:I])/10 .-.05
βstart = rand(m.parLength[:β])/10 .-.05
γstart = rand(m.parLength[:γ]*length(m.data._prodInteract))/10 .-.05
σstart = rand(m.parLength[:σ])/10 .- .05
FEstart = rand(m.parLength[:FE])/100 .-.005

p0 = vcat(Istart,βstart,γstart,σstart,FEstart)
# par0 = parDict(m,p0)
p_est,fval,flag = gradient_ascent(m,p0,grad_tol = 1e-8,max_itr=50)
println("Run NR Estimation")
p_disp = p_est[1:20]
println("Starting at $p_disp")
p_est,fval,flag = newton_raphson_ll(m,p_est)


filename = "flexible_test"
basic_out_3 = MainSpec(df_LA,filename,
    haltonDim = 1,
    spec_per = [:hh_year_id],
    spec_prd = [:product],
    spec_ch = [:choice],
    spec_ch_last = [:iplan],
    spec_prodchr = [:padj,:padj_defpadj],
    spec_prodchr_0= Vector{Symbol}(undef,0),
    spec_inertchr=Vector{Symbol}(undef,0),
    spec_demR=Vector{Symbol}(undef,0),
    spec_fixInt=[:metal,:netname],
    spec_fixEff=[:def_metal,:def_netfe],
    spec_wgt= [:constant],
    method="ga",
    ga_itr=100)

#
# individual_values!(m,par0)
# individual_shares(m,par0)
app = iterate(eachperson(m.data),5)[1]
#
#
ll = log_likelihood(m,par0)

println("Gradient Test")
# p0 = est[1]
par0 = parDict(m,p0)
grad_2 = Vector{Float64}(undef,length(p0))
hess_2 = Matrix{Float64}(undef,length(p0),length(p0))
ll = log_likelihood!(hess_2,grad_2,m,par0)
# # #
# # #
f_obj(x) = log_likelihood(m,x)
grad_1 = Vector{Float64}(undef,length(p0))
hess_1 = Matrix{Float64}(undef,length(p0),length(p0))
fval_old = f_obj(p0)
# # # # # #
println("Grad")
ForwardDiff.gradient!(grad_1,f_obj, p0)#, cfg)
# println("Hessian")
# # cfg = ForwardDiff.HessianConfig(f_obj, p0, ForwardDiff.Chunk{3}())
# # ForwardDiff.hessian!(hess_1,f_obj, p0)
# #
# println(fval_old-ll)
println(maximum(abs.(grad_1-grad_2)))
# println(maximum(abs.(hess_1-hess_2)))

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
Juno.@profile res =  log_likelihood!(hess_2,grad_2,m,par0)
Juno.profiletree()
Juno.profiler()

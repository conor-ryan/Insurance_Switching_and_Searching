import Base.getindex, Base.setindex!, Base.show
# using NLopt
using ForwardDiff


########
# Parameter Structure
#########


mutable struct parDict{T}
    # Parameters
    γ_0::T
    γ::Vector{T}
    I::Vector{T}
    β_0::Vector{T}
    β::Matrix{T}
    σ::Vector{T}
    FE::Matrix{T}
    #Random Coefficients stored (function of σ and draws)
    randCoeffs::Array{T,2}
    # δ values for (ij) pairs
    # δ::Vector{T}
    # Non-delta utility for (ij) pairs and draws
    μ_ij::Matrix{T}
    # Shares for (ij) pairs
    s_hat::Vector{T}
    s_hat_uncond::Vector{T}
    # Search Probability
    ω_i::Vector{T}
    # Total Likelihood
    L_i::Vector{T}

    # r_hat::Vector{T}
    # Share Parameter Derivatives for Products x Parameters
    # dSdθ_j::Matrix{T}
    # dRdθ_j::Matrix{T}

    # d2Sdθ_j::Array{T,3}
    # d2Rdθ_j::Array{T,3}
end

function parDict(m::InsuranceLogit,x::Array{T}) where T
    # Parameter Lengths from model
    #γlen = 1 + m.parLength[:γ]
    γlen = 0
    Ilen = γlen + m.parLength[:I]
    β0len = Ilen + m.parLength[:β]
    βlen = β0len + m.parLength[:γ]*length(m.data._prodInteract)
    σlen = βlen  + m.parLength[:σ]
    FElen = σlen + m.parLength[:FE]

    # Distribute Parameters
    # γ_0 = x[1]
    # γ = x[2:γlen]
    γ_0 = 0.0
    γ = [0.0]
    I_vec = x[(γlen+1):Ilen]
    β_0 = x[(Ilen+1):β0len]
    β_vec = x[(β0len+1):βlen]
    σ_vec = x[(βlen+1):σlen]
    FE_vec = x[(σlen+1):FElen]

    # Store FE as row Vector
    FE = Matrix{T}(undef,1,length(FE_vec))
    FE[1,:] = FE_vec

    # Fill in σ
    σ = Vector{T}(undef,m.parLength[:β])
    σ[:] .= 0.0
    σ[m.data._randCoeffs] = σ_vec

    # Stack Beta into a matrix
    K = m.parLength[:β]
    N = m.parLength[:γ]
    β = Matrix{T}(undef,K,N)
    β[:].=0.0


    ind = 0
    for i in 1:N, j in m.data._prodInteract
        ind+=1
        β[j,i] = β_vec[ind]
    end

    #Calculate Random Coefficients matrix
    (S,R) = size(m.draws)
    randCoeffs = Array{T,2}(undef,S,m.parLength[:β])
    calcRC!(randCoeffs,σ,m.draws,m.data._randCoeffs)

    #Initialize (ij) pairs of deltas
    L, M = size(m.data.data)
    if S>1
        μ_ij = Matrix{T}(undef,S,M)
    else
        μ_ij = Matrix{T}(undef,1,M)
    end
    s_hat = Vector{T}(undef,M)
    s_hat_uncond = Vector{T}(undef,M)

    max_per = Int.(maximum(person(m.data)))
    L_i = Vector{T}(undef,max_per)
    max_srch = maximum(m.data._searchDict[max_per])
    ω_i = Vector{T}(undef,max_srch)


    return parDict{T}(γ_0,γ,I_vec,β_0,β,σ,FE,randCoeffs,μ_ij,s_hat,s_hat_uncond,ω_i,L_i)
end

function calcRC!(randCoeffs::Array{S,2},σ::Array{T,1},draws::Array{Float64,2},randIndex::Vector{Int}) where {T,S}
    (N,K) = size(randCoeffs)
    k_ind = 0
    for k in 1:K
        if k in randIndex
            k_ind += 1
            for n in 1:N
                randCoeffs[n,k] = draws[n,k_ind]*σ[k]
            end
        else
            randCoeffs[:,k] .= 0.0
        end
    end
    return Nothing
end


###########
# Calculating Preferences
###########


function calc_indCoeffs(p::parDict{T},β::Array{T,1},d::S) where {T,S}
    Q = length(β)
    (N,K) = size(p.randCoeffs)
    β_i = Array{T,2}(undef,N,Q)
    γ_i = d
    for n in 1:N
        β_i[n,1] = β[1]
    end

    for k in 1:K, n in 1:N
        β_i[n,k] = β[k] + p.randCoeffs[n,k]
    end

    β_i = permutedims(β_i,(2,1))
    return β_i, γ_i
end

function individual_values!(d::InsuranceLogit,p::parDict{T}) where T
    # Calculate μ_ij, which depends only on parameters
    for app in eachperson(d.data)
        util_value!(app,p)
    end
    return Nothing
end

function util_value!(app::ChoiceData,p::parDict{T}) where T
    γ_0 = p.γ_0
    γ = p.γ
    ι = p.I
    β_0= p.β_0
    M1 = length(β_0)
    β = p.β
    (M2,N2) = size(β)
    fe = p.FE
    randIndex = app._randCoeffs

    ind = person(app)[1]
    idxitr = app._personDict[ind]
    @inbounds X = permutedims(prodchars(app),(2,1))
    @inbounds Z = demoRaw(app)[:,1] # Currently Requiring Demographics are Constant
    X_last_all =inertchars(app)#[:,1]
    y_last = choice_last(app)
    elig = autoelig(app)
    F = fixedEffects(app,idxitr)

    _yearDict = app._personYearDict[ind]
    years = sort(Int.(keys(_yearDict)))
    ret = zeros(length(years))
    X_last = Matrix{Float64}(undef,length(years),size(X_last_all,1))
    for (i,yr) in enumerate(years)
        @inbounds ret[i] = maximum(elig[_yearDict[yr]])
        @inbounds X_last[i,:] = X_last_all[:,_yearDict[yr][1]]
    end


    # FE is a row Vector
    # if T== Float64
    #     controls = zeros(length(idxitr))
    #     for k in 1:length(controls)
    #         for j in app._rel_fe_Dict[ind]
    #             controls[k]+= fe[j]*F[j,k]
    #         end
    #     end
    # else
    controls = fe*F
    # end

    # demos = γ_0 + dot(γ,Z)
    demos = 0.0

    if M2>0
        β_z = β*Z
        β_i, γ_i = calc_indCoeffs(p,β_z,demos)
        chars_int = X*β_i
        # chars = diag(X*β_z)
        # γ_i = 0.0
    else
        chars = zeros(length(idxitr),1)
        γ_i = 0.0
    end

    if M1>0
        chars_0 = X*β_0
    else
        chars_0 = zeros(length(idxitr))
    end

    if length(ι)>0
        search = exp.(X_last*ι)
        s_prob = (1 .- ret) .+ ret.*(search./(1 .+ search))
        p.ω_i[app._searchDict[ind]] = s_prob
    else
        p.ω_i[:].=1.0
    end

    K = length(idxitr)
    N = size(p.randCoeffs,1)
    for k = 1:K,n = 1:N
        @fastmath u = exp(chars_int[k,n] + chars_0[k] + controls[k] + γ_i)
        p.μ_ij[n,idxitr[k]] = u
    end

    return Nothing
end


function calc_shares(μ_ij::Array{T,2},ω_i::Vector{T},iplan::Vector{Float64},y::Vector{Int},dict::Dict{Int,UnitRange}) where T
    (N,K) = size(μ_ij)
    s_hat = Matrix{T}(undef,K,N)
    yr_next = findYearInd(dict)
    expsum = Vector{T}(undef,length(yr_next)+1)
    for n in 1:N
        expsum[:].=0.0
        yr_ind = 1
        for k in 1:K
            if k in yr_next
                yr_ind +=1
            end
            expsum[yr_ind] += μ_ij[n,k]
        end
        yr_ind = 1
        for k in 1:K
            if k in yr_next
                yr_ind +=1
            end
            s_hat[k,n] = iplan[k]*(1-ω_i[yr_ind]) + ω_i[yr_ind]*(μ_ij[n,k]/expsum[yr_ind])
            # s_hat[k,n] =(μ_ij[n,k]/expsum[yr_ind])
        end
    end
    s_mean = mean(s_hat,dims=2)

    chosen = s_hat[y,:]
    ll = prod(chosen,dims=1)
    ll_mean = log(mean(ll))
    # ll_mean = mean(ll)

    return s_mean,ll_mean
end

function calc_shares_mat(μ_ij::Array{T,2},ω_i::Vector{T},iplan::Vector{Float64},y::Vector{Int},yr_next::Vector{Int}) where T
    (N,K) = size(μ_ij)
    s_hat = Matrix{T}(undef,K,N)
    s_cond = Matrix{T}(undef,K,N)
    expsum = Matrix{T}(undef,length(yr_next)+1,N)
    expsum[:].=0.0
    for n in 1:N
        yr_ind = 1
        for k in 1:K
            if k in yr_next
                yr_ind +=1
            end
            expsum[yr_ind,n] += μ_ij[n,k]
        end
        yr_ind = 1
        for k in 1:K
            if k in yr_next
                yr_ind +=1
            end
            s = μ_ij[n,k]/expsum[yr_ind,n]
            s_hat[k,n] = iplan[k]*(1-ω_i[yr_ind]) + ω_i[yr_ind]*s
            s_cond[k,n] = s
        end
    end

    chosen = s_hat[y,:]
    ll = prod(chosen,dims=1)[:]
    # ll_mean = log(mean(ll))

    return s_hat,s_cond,ll,expsum
end

function individual_shares(d::InsuranceLogit,p::parDict{T}) where T
    # Store Parameters
    μ_ij_large = p.μ_ij
    ω_large = p.ω_i
    iplan_large = choice_last(d.data)[:]
    y_large = choice(d.data)[:]

    for (ind,idxitr) in d.data._personDict
        # println(ind)
        u = μ_ij_large[:,idxitr]
        ω = ω_large[d.data._searchDict[ind]]
        iplan = iplan_large[idxitr]
        dict = d.data._personYearDict[ind]
        y = y_large[idxitr]
        y = findall(y.==1)
        s_hat,ll = calc_shares(u,ω,iplan,y,dict)
        p.s_hat[idxitr] = s_hat
        p.L_i[ind] = ll
    end
    return Nothing
end


function ind_wtp(app::ChoiceData,p::parDict{T}) where T
    γ_0 = p.γ_0
    γ = p.γ
    ι = p.I
    β_0= p.β_0
    M1 = length(β_0)
    β = p.β
    (M2,N2) = size(β)
    fe = p.FE
    randIndex = app._randCoeffs

    ind = person(app)[1]
    wgt = weight(app)[1]
    idxitr = app._personDict[ind]
    X = permutedims(prodchars(app),(2,1))
    Z = demoRaw(app)[:,1]
    X_last =inertchars(app)[:,1]
    y_last = choice_last(app)



    # demos = γ_0 + dot(γ,Z)
    demos = 0.0

    β_z = β*Z
    β_i, γ_i = calc_indCoeffs(p,β_z,demos)

    β_mean = mean(β_i,dims=2) + β_0
    return β_mean, wgt
end


function coeff_values(d::InsuranceLogit,p::parDict{T}) where T
    # Calculate μ_ij, which depends only on parameters
    coeff_mat = Matrix{Float64}(undef,length(d.data._personIDs),d.parLength[:β])
    i = 0
    for app in eachperson(d.data)
        i+=1
        β_mean,wgt = ind_wtp(app,p)
        coeff_mat[i,:] = β_mean[:]
    end
    return coeff_mat
end

function marginalEffects(d::InsuranceLogit,p::parDict{T}) where T
    # Calculate μ_ij, which depends only on parameters
    dω_i = Matrix{Float64}(undef,length(p.ω_i),d.parLength[:I])
    returning = Vector{Float64}(undef,length(p.ω_i))
    for app in eachperson(d.data)
        indMargEffect(dω_i,returning,app,p)
    end
    dω_i[:,6] = 10*p.I[6].*p.ω_i.*(1 .- p.ω_i)
    dω_i = dω_i[findall(returning.==1),:]
    ### Manually do marginal effect for dprem
    return round.(100 .*mean(dω_i,dims=1)[2:6],digits=2)
end

function indMargEffect(dω_i::Matrix{Float64},returning::Vector{Float64},app::ChoiceData,p::parDict{T}) where T
    γ_0 = p.γ_0
    γ = p.γ
    ι = p.I
    β_0= p.β_0
    M1 = length(β_0)
    β = p.β
    (M2,N2) = size(β)
    fe = p.FE
    randIndex = app._randCoeffs

    ind = person(app)[1]
    idxitr = app._personDict[ind]
    @inbounds X = permutedims(prodchars(app),(2,1))
    @inbounds Z = demoRaw(app)[:,1] # Currently Requiring Demographics are Constant
    X_last_all =inertchars(app)#[:,1]
    y_last = choice_last(app)
    elig = autoelig(app)
    F = fixedEffects(app,idxitr)

    _yearDict = app._personYearDict[ind]
    years = sort(Int.(keys(_yearDict)))
    ret = zeros(length(years))
    X_last = Matrix{Float64}(undef,length(years),size(X_last_all,1))
    for (i,yr) in enumerate(years)
        @inbounds ret[i] = maximum(elig[_yearDict[yr]])
        @inbounds X_last[i,:] = X_last_all[:,_yearDict[yr][1]]
    end

    (N,K) = size(X_last)

    for k in 1:K
        X_last0 = copy(X_last)
        X_last0[:,k] .= 0.0
        X_last1 = copy(X_last)
        X_last1[:,k] .= 1.0
        search0 = exp.(X_last0*ι)
        search1 = exp.(X_last1*ι)
        search_prob0 = search0./(1 .+ search0)
        search_prob1 = search1./(1 .+ search1)
        dω_i[app._searchDict[ind],k] = ret.*(search_prob1-search_prob0)
    end
    returning[app._searchDict[ind]] = ret[:]
    return nothing
end

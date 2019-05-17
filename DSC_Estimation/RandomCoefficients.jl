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

    ω_i = Vector{T}(undef,Int.(maximum(person(m.data))))

    # Deltas are turned off
    # δ = Vector{T}(undef,M)
    # unpack_δ!(δ,m)
    # δ = ones(M)

    # Q = m.parLength[:All]
    # J = length(m.prods)
    # dSdθ_j = Matrix{T}(undef,Q,J)
    # dRdθ_j = Matrix{T}(undef,Q,J)
    # d2Sdθ_j = Array{T,3}(undef,Q,Q,J)
    # d2Rdθ_j = Array{T,3}(undef,Q,Q,J)

    return parDict{T}(γ_0,γ,I_vec,β_0,β,σ,FE,randCoeffs,μ_ij,s_hat,s_hat_uncond,ω_i)
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

# function calc_indCoeffs{T}(p::parDict{T},β::Array{T,1},d::T)
#     Q = length(β)
#     (N,K) = size(p.randCoeffs)
#     β_i = Array{T,2}(N,Q)
#     γ_i = d
#     β_i[:,1] = β[1]
#
#     for k in 2:Q, n in 1:N
#         β_i[n,k] = β[k] + p.randCoeffs[n,k-1]
#     end
#
#     β_i = permutedims(β_i,(2,1))
#     return β_i, γ_i
# end

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
    X = permutedims(prodchars(app),(2,1))
    Z = demoRaw(app)[:,1]
    X_last =inertchars(app)[:,1]
    y_last = choice_last(app)
    F = fixedEffects(app,idxitr)

    # FE is a row Vector
    if T== Float64
        controls = zeros(length(idxitr))
        for k in 1:length(controls)
            for j in app._rel_fe_Dict[ind]
                controls[k]+= fe[j]*F[j,k]
            end
        end
    else
        controls = fe*F
    end

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

    search = exp(dot(X_last,ι))
    if (length(ι)>0) & (sum(y_last)>0)
        p.ω_i[Int.(ind)] = search/(1+search)
    else
        p.ω_i[Int.(ind)] = 1.0
    end


    K = length(idxitr)
    N = size(p.randCoeffs,1)
    for k = 1:K,n = 1:N
        @fastmath u = exp(chars_int[k,n] + chars_0[k] + controls[k] + γ_i)
        p.μ_ij[n,idxitr[k]] = u
    end

    return Nothing
end


function calc_shares(μ_ij::Array{T},ω_i::T,iplan::Array{Float64}) where T
    (N,K) = size(μ_ij)
    s_hat = Matrix{T}(undef,K,N)
    for n in 1:N
        expsum = 0.0
        for i in 1:K
            expsum += μ_ij[n,i]
        end
        for i in 1:K
            s_hat[i,n] = μ_ij[n,i]/expsum
        end
    end
    s_mean_uncond = mean(s_hat,dims=2)
    s_hat = similar(s_mean_uncond)
    for i in 1:K
        s_hat[i] = iplan[i]*(1-ω_i) + ω_i*(s_mean_uncond[i])
    end
    return s_hat,s_mean_uncond
end

function individual_shares(d::InsuranceLogit,p::parDict{T}) where T
    # Store Parameters
    μ_ij_large = p.μ_ij
    ω_large = p.ω_i
    iplan_large = choice_last(d.data)[:]
    for (ind,idxitr) in d.data._personDict
        u = μ_ij_large[:,idxitr]
        ω = ω_large[ind]
        iplan = iplan_large[idxitr]
        s_hat,s_uncond = calc_shares(u,ω,iplan)
        p.s_hat[idxitr] = s_hat
        p.s_hat_uncond[idxitr] = s_uncond
    end
    return Nothing
end

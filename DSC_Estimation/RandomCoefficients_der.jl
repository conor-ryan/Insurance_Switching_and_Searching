function unPackChars(app::ChoiceData,d::InsuranceLogit)
    ind = person(app)[1]
    S_ij = transpose(choice(app))
    wgt = transpose(weight(app))
    idxitr = d.data._personDict[ind]

    X_t = prodchars(app)
    X_0_t = prodchars0(app)
    X_int = X_t[d.data._prodInteract,:]
    X_last = inertchars(app)
    y_last = choice_last(app)[:]
    Z = demoRaw(app)[:,1]
    F_t = fixedEffects(app,idxitr)

    return ind, S_ij, wgt,idxitr, X_t, X_0_t,X_int, Z, F_t, X_last, y_last
end

function unPackParChars(p::parDict{T},idxitr::UnitRange{Int}) where T
    μ_ij = p.μ_ij[:,idxitr]
    #Get Market Shares
    s_hat = p.s_hat[idxitr]
    s_uncond = p.s_hat_uncond[idxitr]
    return μ_ij,s_hat,s_uncond
end

function sumShares!(s_hat::Vector{T}) where T
    # Fix possible computational error
    for k in eachindex(s_hat)
        if abs(s_hat[k])<=1e-300
            s_hat[k]=1e-15
            #println("Hit Share Constraint for person $ind, product $k")
        end
    end
    # s_insured = sum(s_hat)
    # if s_insured>=(1-1e-300)
    #     s_insured= 1 - 1e-15
    #     #println("Hit insured constraint for person $ind")
    # end

    return nothing
end

function relPar(app::ChoiceData,d::InsuranceLogit,F_t::SubArray,ind::Float64)
    Q = d.parLength[:All]
    Q_0 = Q - size(F_t,1)
    ## Relevant Parameters for this observation
    pars_relevant = vcat(1:Q_0,Q_0 .+ app._rel_fe_Dict[ind])

    return pars_relevant
end


function preCalcμ(μ_ij::Matrix{T},dict::Dict{Int,UnitRange}) where T
    (N,K) = size(μ_ij)
    yr_next = findYearInd(dict)
    expsum = zeros(length(yr_next)+1,N)
    for n in 1:N
        yr_ind = 1
        for k in 1:K
            if k in yr_next
                yr_ind +=1
            end
            expsum[yr_ind,n] += μ_ij[n,k]
        end
    end
    return expsum
end

function findYearInd(dict)
    years = sort(Int.(keys(dict)))
    if length(years)>1
        yr_next = Vector{Int}(undef,length(years)-1)
        for (i,yr) in enumerate(years[2:end])
            yr_next[i] = minimum(dict[yr])
        end
    else
        yr_next = Vector{Int}(undef,0)
    end
    return yr_next
end

function ll_Terms(wgt::Array{Float64,N},S_ij::Array{Float64,N},
                    s_hat::Vector{T}) where {N,T}
    K = length(S_ij)
    ll_obs = 0.0
    gll_t1 = Vector{T}(undef,K)
    gll_t2 = Vector{T}(undef,K)
    gll_t3 = Vector{T}(undef,K)
    # gll_t4 = Vector{T}(undef,K)
    # gll_t5 = Vector{T}(undef,K)
    # gll_t6 = Vector{T}(undef,K)
    for k in 1:K
        #Gradient Terms
        gll_t1[k] = wgt[k]*S_ij[k]*(1/s_hat[k])
        gll_t2[k] = wgt[k]*S_ij[k]*(1/s_hat[k]^2)
        gll_t3[k] = wgt[k]*S_ij[k]*2*(1/s_hat[k]^3)
        # gll_t4[k] = wgt[k]*S_ij[k]*urate[k]*(1/(s_insured) + 1/(1-s_insured))
        # gll_t5[k] = wgt[k]*S_ij[k]*urate[k]*(1/(s_insured^2) - 1/((1-s_insured)^2))
        # gll_t6[k] = wgt[k]*S_ij[k]*urate[k]*2*(1/(s_insured^3) + 1/((1-s_insured)^3))
        # Log Likelihood
        # ll_obs+=wgt[k]*S_ij[k]*(log(s_hat[k]) -
        #                 urate[k]*(log(s_insured)-log(1-s_insured)))
        ll_obs+=wgt[k]*S_ij[k]*log(s_hat[k])
    end

    return gll_t1, gll_t2, gll_t3, ll_obs
end



function ll_obs!(hess::Matrix{Float64},grad::Vector{Float64},
                    app::ChoiceData,d::InsuranceLogit,p::parDict{T}) where T

        ind, S_ij, wgt, idxitr, X_t, X_0_t,X_int, Z, F_t, X_last, y_last = unPackChars(app,d)
        wgt = convert(Array{Float64,2},wgt)
        S_ij = convert(Array{Float64,2},S_ij)
        ω_i = p.ω_i[app._searchDict[ind]]
        choice_ind = findall(S_ij[:].==1)

        prodidx = Int.(product(app))

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat,s_uncond = unPackParChars(p,idxitr)

        sumShares!(s_hat)

        (N,K) = size(μ_ij)

        # Initialize Gradient
        #(Q,N,K) = size(dμ_ij)
        pars_relevant = d._rel_par_Dict[ind]
        # pars_relevant = 1:d.parLength[:All]

        # Year Indices
        yr_next = findYearInd(app._personYearDict[ind])

        # Pre-Calculate Shares
        # μ_ij_sums = preCalcμ(μ_ij,app._personYearDict[ind])
        s_hat, s_n, ll_n, μ_ij_sums = calc_shares_mat(μ_ij,ω_i,y_last,choice_ind,yr_next)
        ω_yr = expandYear(ω_i,yr_next,K)
        ll_mean = log(mean(ll_n))

        # Pre-Calculate Log-Likelihood Terms for Gradient
        # Also Calculate Log-Likelihood itself
        # gll_t1, gll_t2, gll_t3, ll_obs = ll_Terms(wgt,S_ij,s_hat)

        #hess = zeros(Q,Q)
        #hess[:] = 0.0
        #grad[:] = 0.0
        X_mat = Array{Float64}(undef,N,K)
        Y_list = Array{Union{Float64, Array{Float64,1}, Array{Float64,2}},1}(undef,length(pars_relevant))
        #Y_mat = Array{Float64}(N,K)

        # Allocate Memory

        dS_x = Matrix{Float64}(undef,length(choice_ind),N)
        dS_x_full = Matrix{Float64}(undef,length(choice_ind),N)
        dS_xy = Matrix{Float64}(undef,length(choice_ind),N)
        dlnL_list = Vector{Vector{Float64}}(undef,length(pars_relevant))


        dS_x_list = Vector{Matrix{Float64}}(undef,length(pars_relevant))
        dS_full_list = Vector{Matrix{Float64}}(undef,length(pars_relevant))

        #γlen = 1 + d.parLength[:γ]
        γlen = 0
        Ilen = γlen + d.parLength[:I]
        β0len = Ilen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]*length(d.data._prodInteract)
        σlen = βlen  + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]

        s_adjust = s_uncond .- y_last

        ll = 0
        for (q_i,q) in enumerate(pars_relevant)
            X = returnParameterX!(q,X_mat,
                            Z,X_0_t,X_t,X_last,X_int,draws,F_t,
                            γlen,Ilen,β0len,βlen,σlen)
            @inbounds Y_list[q_i] = copy(X)

            if (γlen < q <= Ilen)
                dlldβ, ll, dlnL_x = grad_calc_ω!(dS_x,
                            s_n,s_hat,ll_n,
                            X,
                            yr_next,
                            choice_ind,
                            y_last,
                            ω_yr)
                (K,N) = size(dS_x)
                for (c,k) in enumerate(choice_ind)
                    for n in 1:N
                    dS_x_full[c,n] = dS_x[c,n]*(s_n[k,n] - y_last[k])
                    end
                end
            else
                dlldβ, ll, dlnL_x = grad_calc_s!(dS_x,
                            s_n,s_hat,ll_n,
                            μ_ij,X,
                            μ_ij_sums,
                            yr_next,
                            choice_ind,
                            ω_yr)
                (K,N) = size(dS_x)
                for k in 1:K, n in 1:N
                    dS_x_full[k,n] = dS_x[k,n]*ω_i[k]
                end
            end
            # Save Answers
            dS_x_list[q_i] = dS_x[:,:]
            dS_full_list[q_i] = dS_x_full[:,:]
            dlnL_list[q_i] = dlnL_x[:]

            ## Calculate Gradient
            grad[q]+= dlldβ


            for (r_i,r) in enumerate(pars_relevant)
                if r>q
                    continue
                end
                Y_mat = Y_list[r_i]


                dS_y = dS_x_list[r_i]
                dS_y_full = dS_full_list[r_i]
                dlnL_y = dlnL_list[r_i]
                # dω_y = dω_list[r_i]
                # dS_u_y = dS_uncond_list[r_i]

                if (q <= Ilen)
                    hess_calc_ω!(dS_xy,
                                s_n,
                                X,Y_mat,
                                yr_next,
                                choice_ind,
                                y_last,
                                ω_yr)
                elseif (r <= Ilen)
                    dS_xy = dS_y.*dS_x
                else
                    hess_calc_s!(dS_xy,
                                s_n,
                                μ_ij,X,Y_mat,
                                μ_ij_sums,
                                yr_next,
                                choice_ind,
                                ω_yr)
                end

                hess_obs = hess_calc_ll(dS_xy,dS_x_full,dS_y_full,
                                        dlnL_x,dlnL_y,ll_n,
                                        s_hat,choice_ind)

                hess[q,r]+= hess_obs
                if (q!=r)
                    hess[r,q]+= hess_obs
                end
            end
        end

    return ll_mean
end


function ll_obs!(grad::Vector{Float64},
                    app::ChoiceData,d::InsuranceLogit,p::parDict{T}) where T

        ind, S_ij, wgt, idxitr, X_t, X_0_t,X_int, Z, F_t, X_last, y_last = unPackChars(app,d)
        wgt = convert(Array{Float64,2},wgt)
        S_ij = convert(Array{Float64,2},S_ij)
        ω_i = p.ω_i[app._searchDict[ind]]
        choice_ind = findall(S_ij[:].==1)

        prodidx = Int.(product(app))

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat,s_uncond = unPackParChars(p,idxitr)

        sumShares!(s_hat)

        (N,K) = size(μ_ij)

        # Initialize Gradient
        #(Q,N,K) = size(dμ_ij)
        pars_relevant = d._rel_par_Dict[ind]
        # pars_relevant = 1:d.parLength[:All]

        # Year Indices
        yr_next = findYearInd(app._personYearDict[ind])

        # Pre-Calculate Shares
        # μ_ij_sums = preCalcμ(μ_ij,app._personYearDict[ind])
        s_hat, s_n, ll_n, μ_ij_sums = calc_shares_mat(μ_ij,ω_i,y_last,choice_ind,yr_next)
        ω_yr = expandYear(ω_i,yr_next,K)
        ll_mean = log(mean(ll_n))

        # Pre-Calculate Log-Likelihood Terms for Gradient
        # Also Calculate Log-Likelihood itself
        # gll_t1, gll_t2, gll_t3, ll_obs = ll_Terms(wgt,S_ij,s_hat)

        #hess = zeros(Q,Q)
        #hess[:] = 0.0
        #grad[:] = 0.0
        X_mat = Array{Float64}(undef,N,K)
        # Y_list = Array{Union{Float64, Array{Float64,1}, Array{Float64,2}},1}(undef,length(pars_relevant))
        #Y_mat = Array{Float64}(N,K)

        # Allocate Memory

        dS_x = Matrix{Float64}(undef,length(choice_ind),N)
        # dS_x_full = Matrix{Float64}(undef,length(choice_ind),N)
        # dS_xy = Matrix{Float64}(undef,length(choice_ind),N)
        # dlnL_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        #
        #
        # dS_x_list = Vector{Matrix{Float64}}(undef,length(pars_relevant))
        # dS_full_list = Vector{Matrix{Float64}}(undef,length(pars_relevant))

        #γlen = 1 + d.parLength[:γ]
        γlen = 0
        Ilen = γlen + d.parLength[:I]
        β0len = Ilen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]*length(d.data._prodInteract)
        σlen = βlen  + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]

        s_adjust = s_uncond .- y_last

        ll = 0
        for (q_i,q) in enumerate(pars_relevant)
            X = returnParameterX!(q,X_mat,
                            Z,X_0_t,X_t,X_last,X_int,draws,F_t,
                            γlen,Ilen,β0len,βlen,σlen)
            # @inbounds Y_list[q_i] = copy(X)

            if (γlen < q <= Ilen)
                dlldβ, ll, dlnL_x = grad_calc_ω!(dS_x,
                            s_n,s_hat,ll_n,
                            X,
                            yr_next,
                            choice_ind,
                            y_last,
                            ω_yr)
            else
                dlldβ, ll, dlnL_x = grad_calc_s!(dS_x,
                            s_n,s_hat,ll_n,
                            μ_ij,X,
                            μ_ij_sums,
                            yr_next,
                            choice_ind,
                            ω_yr)
            end

            ## Calculate Gradient
            grad[q]+= dlldβ
        end

    return ll_mean
end

function hess_calc_ll(dS_xy::Matrix{T},dS_x::Matrix{T},dS_y::Matrix{T},
                    dlnLd_x::Vector{T},dlnLd_y::Vector{T},ll_n::Vector{T},
                    s_hat::Matrix{T},choice_ind::Vector{Int}) where {T}

    N = length(ll_n)
    LL = 0.0
    d2Ldβ = 0.0
    dLdβ_x = 0.0
    dLdβ_y = 0.0
    for n in 1:N
        d2lnLdβ = calc_d2ll(n,dS_xy,dS_x,dS_y,s_hat,choice_ind)
        LL += ll_n[n]
        d2Ldβ += ll_n[n]*(d2lnLdβ+dlnLd_x[n]*dlnLd_y[n])
        dLdβ_x += dlnLd_x[n]*ll_n[n]
        dLdβ_y += dlnLd_y[n]*ll_n[n]
    end
    dLLdx = dLdβ_x/LL
    dLLdy = dLdβ_y/LL
    d2LLdβ = d2Ldβ/LL - dLLdx*dLLdy

    return d2LLdβ
end


function grad_calc_s!(dS_x::Matrix{T},s_n::Matrix{T},s_hat::Matrix{T},
                    ll_n::Vector{Float64},
                    μ_ij::Matrix{T},
                    X_mat::S,
                    μ_ij_sums::Matrix{T},
                    yr_next::Vector{Int},
                    choice_ind::Vector{Int},
                    ω_i::Vector{Float64}) where {S,T}

    (N,K) = size(μ_ij)
    dμ_ij_x_sums = Vector{T}(undef,length(yr_next)+1)

    dLdβ = 0.0
    LL = 0.0
    Γ_x = Vector{T}(undef,length(yr_next)+1)
    dlnLdβ = Vector{T}(undef,N)
    for n in 1:N
        # μ_ij_sums_n = μ_ij_sums[:,n]
        calc_derSums_x!(n,dμ_ij_x_sums,X_mat,μ_ij,yr_next)

        for y in 1:length(Γ_x)
            Γ_x[y] = dμ_ij_x_sums[y]/μ_ij_sums[y,n]
        end
        # Γ_x = expandYear(Γ_x,yr_next,K)

        # s_hat = y_last.*(1 .- ω_yr) + ω_yr.*s_n

        dl = calc_dll(dS_x,n,X_mat,Γ_x,ω_i,s_n,s_hat,choice_ind)
        dlnLdβ[n] = dl
        LL += ll_n[n]
        dLdβ += dl*ll_n[n]
        # s_n = s_n./μ_ij_sums_n
        # calc_prodTerms_x!(n,dS_x,X_mat,s_n,Γ_x)
    end

    dLLdβ = dLdβ/LL
    LL = log(LL/N)
    return dLLdβ, LL, dlnLdβ
end



function hess_calc_s!(dS_xy::Matrix{T},s_n::Matrix{T},
                    μ_ij::Matrix{T},
                    X_mat::S,Y_mat::R,
                    μ_ij_sums::Matrix{T},
                    yr_next::Vector{Int},
                    choice_ind::Vector{Int},
                    ω_i::Vector{Float64}) where {S,R,T}

    (N,K) = size(μ_ij)
    dμ_ij_x_sums = Vector{T}(undef,length(yr_next)+1)
    dμ_ij_y_sums = Vector{T}(undef,length(yr_next)+1)
    dμ_ij_xy_sums = Vector{T}(undef,length(yr_next)+1)

    Γ_x = Vector{T}(undef,length(yr_next)+1)
    Γ_y = Vector{T}(undef,length(yr_next)+1)
    Γ_xy = Vector{T}(undef,length(yr_next)+1)
    for n in 1:N
        calc_derSums_xy!(n,dμ_ij_x_sums,dμ_ij_y_sums,dμ_ij_xy_sums,
                        X_mat,Y_mat,μ_ij,yr_next)

        for y in 1:length(Γ_x)
            Γ_x[y] = dμ_ij_x_sums[y]/μ_ij_sums[y,n]
            Γ_y[y] = dμ_ij_y_sums[y]/μ_ij_sums[y,n]
            Γ_xy[y] = dμ_ij_xy_sums[y]/μ_ij_sums[y,n]
        end

        calc_d2S(dS_xy,n,X_mat,Y_mat,Γ_x,Γ_y,Γ_xy,ω_i,s_n,choice_ind)
    end
    return nothing
end

function grad_calc_ω!(dS_x::Matrix{T},s_n::Matrix{T},s_hat::Matrix{T},
                    ll_n::Vector{Float64},
                    X_mat::S,
                    yr_next::Vector{Int},
                    choice_ind::Vector{Int},
                    y_last::Vector{Float64},
                    ω_i::Vector{Float64}) where {S,T}

    (K,N) = size(s_n)
    dLdβ = 0.0
    LL = 0.0
    dlnLdβ = Vector{T}(undef,N)
    for n in 1:N
        @inbounds ll = ll_n[n]
        dl = calc_dll_ω(dS_x,n,X_mat,y_last,ω_i,s_n,s_hat,choice_ind)
        dlnLdβ[n] = dl
        LL += ll_n[n]
        dLdβ += dl*ll_n[n]
    end

    dLLdβ = dLdβ/LL
    LL = log(LL/N)
    return dLLdβ, LL, dlnLdβ
end


function hess_calc_ω!(dS_xy::Matrix{T},s_n::Matrix{T},
                    X_mat::S,
                    Y_mat::S,
                    yr_next::Vector{Int},
                    choice_ind::Vector{Int},
                    y_last::Vector{Float64},
                    ω_i::Vector{Float64}) where {S,T}

    (K,N) = size(s_n)
    for n in 1:N
        calc_d2S_ω(dS_xy,n,X_mat,Y_mat,y_last,ω_i,s_n,choice_ind)
    end
    return nothing
end


function expandYear(x::Vector{T},yr_next::Vector{Int},last::Int) where T
    y = Vector{T}(undef,last)
    init = 1
    for (i,n) in enumerate(yr_next)
        y[init:(n-1)] .= x[i]
        init = n
    end
    y[init:last].=x[end]
    return y
end




function returnParameter!(q::Int64,X_mat::Matrix{Float64},
                        Z::Vector{Float64},X_0_t::Matrix{Float64},
                        X_t::Matrix{Float64},X_last::Matrix{Float64},
                        X_int::Matrix{Float64},
                        draws::Matrix{Float64},
                        F_t::SubArray{Float64,2},
                        γlen::Int64,Ilen::Int64,β0len::Int64,βlen::Int64,σlen::Int64)
    (N,K) = size(X_mat)
    (Q,R) = size(X_int)
    if q<0
        X_mat[:] .= 1.0
    elseif q<=γlen
        X_mat[:] .= Z[q]
    elseif q<=Ilen
        for n in 1:N
            @inbounds X_mat[n,:] = X_last[q-γlen,:]
        end
    elseif q<=β0len
        for n in 1:N
            @inbounds X_mat[n,:] = X_t[q-Ilen,:]
        end
    elseif q<=βlen
        # Characteristic Interactions
        Z_ind = Int(ceil((q-β0len)/Q))
        # println(Z_ind)
        X_ind = ((q-β0len-1)%Q)+1
        # println(X_ind)
        for n in 1:N
            @inbounds X_mat[n,:] = X_int[X_ind,:].*Z[Z_ind]
        end
    elseif q<=σlen
        #Quality Random Effect
        for n in 1:N,k in 1:K
            @inbounds X_mat[n,k] = draws[n,q-(βlen)]*X_0_t[q-(βlen),k]
        end
    else
        #Fixed Effect
        # println(q)
        for n in 1:N
            @inbounds X_mat[n,:] = F_t[q-σlen,:]
        end
    end
    return Nothing
end


function returnParameterX!(q::Int64,X_mat::Matrix{Float64},
                        Z::Vector{Float64},X_0_t::Matrix{Float64},
                        X_t::Matrix{Float64},X_last::Matrix{Float64},
                        X_int::Matrix{Float64},
                        draws::Matrix{Float64},
                        F_t::SubArray{Float64,2},
                        γlen::Int64,Ilen::Int64,β0len::Int64,βlen::Int64,σlen::Int64)
    (N,K) = size(X_mat)
    (Q,R) = size(X_int)
    if q<0
        X = 1.0
    elseif q<=γlen
        X = Z[q]
    elseif q<=Ilen
        X = X_last[q-γlen,:]
    elseif q<=β0len
        X = X_t[q-Ilen,:]
    elseif q<=βlen
        # Characteristic Interactions
        Z_ind = Int(ceil((q-β0len)/Q))
        # println(Z_ind)
        X_ind = ((q-β0len-1)%Q)+1
        # println(X_ind)
        X = X_int[X_ind,:].*Z[Z_ind]
    elseif q<=σlen
        #Quality Random Effect
        for n in 1:N,k in 1:K
            @inbounds X_mat[n,k] = draws[n,q-(βlen)]*X_0_t[q-(βlen),k]
        end
        X = copy(X_mat)
    else
        X = F_t[q-σlen,:]
    end
    return X
end


function find_par_rel(app::ChoiceData,d::InsuranceLogit)

        ind, S_ij, wgt, idxitr, X_t, X_0_t,X_int, Z, F_t, X_last, y_last = unPackChars(app,d)
        draws = d.draws

        Q = d.parLength[:All]
        N = size(d.draws,1)
        K = size(X_0_t,2)
        X_mat = Array{Float64}(undef,N,K)

        γlen = 0
        Ilen = γlen + d.parLength[:I]
        β0len = Ilen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]*length(d.data._prodInteract)
        σlen = βlen  + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]

        pars_rel = Int.(1:Q)
        pars_bool = Vector{Bool}(undef,Q)
        ll = 0
        for q in 1:Q
            X = returnParameterX!(q,X_mat,
                            Z,X_0_t,X_t,X_last,X_int,draws,F_t,
                            γlen,Ilen,β0len,βlen,σlen)
            rel = any(abs.(X).>0)
            pars_bool[q] = rel
        end
    return pars_rel[pars_bool]
end

function calc_relParDict!(d::InsuranceLogit)
    for app in eachperson(d.data)
        ind = person(app)[1]
        d._rel_par_Dict[ind] = find_par_rel(app,d)
    end
    return nothing
end

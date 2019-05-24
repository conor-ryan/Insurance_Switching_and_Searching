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


function preCalcμ(μ_ij::Matrix{T}) where T
    μ_ij_sums = sum(μ_ij,dims=2)
    return μ_ij_sums[:]
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
        ω_i = p.ω_i[Int.(ind)]

        prodidx = Int.(product(app))

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat,s_uncond = unPackParChars(p,idxitr)

        sumShares!(s_hat)

        (N,K) = size(μ_ij)

        # Initialize Gradient
        #(Q,N,K) = size(dμ_ij)
        pars_relevant = relPar(app,d,F_t,ind)

        # Pre-Calculate Squares
        μ_ij_sums = preCalcμ(μ_ij)

        # Pre-Calculate Log-Likelihood Terms for Gradient
        # Also Calculate Log-Likelihood itself
        gll_t1, gll_t2, gll_t3, ll_obs = ll_Terms(wgt,S_ij,s_hat)

        #hess = zeros(Q,Q)
        #hess[:] = 0.0
        #grad[:] = 0.0
        X_mat = Array{Float64}(undef,N,K)
        Y_list = Array{Array{Float64,2},1}(undef,length(pars_relevant))
        #Y_mat = Array{Float64}(N,K)

        # Allocate Memory
        # dS_xyz = Vector{Float64}(undef,K)
        dS_xy = Vector{Float64}(undef,K)
        dS_u_x = Vector{Float64}(undef,K)
        dS_x = Vector{Float64}(undef,K)
        dω_x = Vector{Float64}(undef,K)


        dS_x_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        dS_uncond_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        dω_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        # dS_xy_list = Array{Vector{Float64},2}(undef,length(pars_relevant),length(pars_relevant))

        s_n = Vector{Float64}(undef,K)


        #γlen = 1 + d.parLength[:γ]
        γlen = 0
        Ilen = γlen + d.parLength[:I]
        β0len = Ilen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]*length(d.data._prodInteract)
        σlen = βlen  + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]

        s_adjust = s_uncond .- y_last


        for (q_i,q) in enumerate(pars_relevant)
            returnParameter!(q,X_mat,
                            Z,X_0_t,X_t,X_last,X_int,draws,F_t,
                            γlen,Ilen,β0len,βlen,σlen)
            @inbounds Y_list[q_i] = X_mat[:,:]

            if (γlen < q <= Ilen)
                grad_calc_ω!(dω_x,ω_i,X_mat)


                dω_list[q_i] = dω_x[:]
                dS_u_x[:] .= 0.0
                dS_uncond_list[q_i] = dS_u_x[:]

                # @inbounds @fastmath @simd for k in 1:K
                #     dS_x[k] = dω_x[k]*(s_uncond[k] - y_last[k])
                # end
                dS_x[:] = dω_x.*s_adjust
            else
                grad_calc_s!(dS_u_x,s_n,
                            μ_ij,X_mat,
                            μ_ij_sums)

                dS_uncond_list[q_i] = dS_u_x[:]
                dω_x[:] .= 0.0
                dω_list[q_i] = dω_x[:]

                dS_x[:] = (dS_u_x[:].*ω_i)./N
            end
            # Save Answers
            dS_x_list[q_i] = dS_x[:]

            ## Calculate Gradient
            grad[q]+= combine_grad(gll_t1,dS_x)


            for (r_i,r) in enumerate(pars_relevant)
                if r>q
                    continue
                end
                Y_mat = Y_list[r_i]


                dS_y = dS_x_list[r_i]
                dω_y = dω_list[r_i]
                dS_u_y = dS_uncond_list[r_i]

                if (q <= Ilen)
                    hess_calc_ω!(dS_xy,s_uncond,ω_i,
                                X_mat,Y_mat,y_last)
                elseif (r <= Ilen)
                    for k in 1:K
                        dS_xy[k] = ((dS_u_x[k])*(dω_y[k]))./N
                    end
                    # dS_xy[:] = (dω_list[q_i]+dS_uncond_list[q_i]).*(dω_list[r_i]+dS_uncond_list[r_i])
                else
                    hess_calc!(dS_xy,s_n,
                                μ_ij,X_mat,Y_mat,
                                μ_ij_sums)
                    for k in 1:K
                        dS_xy[k] = (dS_xy[k]*ω_i)./N
                    end
                    # dS_xy[:].*=ω_i
                end
                # dS_xy_list[q_i,r_i] = dS_xy[:]

                hess_obs = combine_hess(gll_t1,gll_t2,
                            dS_xy,dS_x,dS_y)

                hess[q,r]+= hess_obs
                if (q!=r)
                    hess[r,q]+= hess_obs
                end
            end
        end

    return ll_obs,pars_relevant
end

function ll_obs!(grad::Vector{Float64},
                    app::ChoiceData,d::InsuranceLogit,p::parDict{T}) where T

        ind, S_ij, wgt, idxitr, X_t, X_0_t,X_int, Z, F_t, X_last, y_last = unPackChars(app,d)
        wgt = convert(Array{Float64,2},wgt)
        S_ij = convert(Array{Float64,2},S_ij)
        ω_i = p.ω_i[Int.(ind)]

        prodidx = Int.(product(app))

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat,s_uncond = unPackParChars(p,idxitr)

        sumShares!(s_hat)

        (N,K) = size(μ_ij)

        # Initialize Gradient
        #(Q,N,K) = size(dμ_ij)
        pars_relevant = relPar(app,d,F_t,ind)

        # Pre-Calculate Squares
        μ_ij_sums = preCalcμ(μ_ij)

        # Pre-Calculate Log-Likelihood Terms for Gradient
        # Also Calculate Log-Likelihood itself
        gll_t1, gll_t2, gll_t3, ll_obs = ll_Terms(wgt,S_ij,s_hat)

        #hess = zeros(Q,Q)
        #hess[:] = 0.0
        #grad[:] = 0.0
        X_mat = Array{Float64}(undef,N,K)
        Y_list = Array{Array{Float64,2},1}(undef,length(pars_relevant))
        #Y_mat = Array{Float64}(N,K)

        # Allocate Memory
        # dS_xyz = Vector{Float64}(undef,K)
        dS_xy = Vector{Float64}(undef,K)
        dS_u_x = Vector{Float64}(undef,K)
        dS_x = Vector{Float64}(undef,K)
        dω_x = Vector{Float64}(undef,K)


        dS_x_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        dS_uncond_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        dω_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        # dS_xy_list = Array{Vector{Float64},2}(undef,length(pars_relevant),length(pars_relevant))

        s_n = Vector{Float64}(undef,K)


        #γlen = 1 + d.parLength[:γ]
        γlen = 0
        Ilen = γlen + d.parLength[:I]
        β0len = Ilen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]*length(d.data._prodInteract)
        σlen = βlen  + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]

        s_adjust = s_uncond .- y_last


        for (q_i,q) in enumerate(pars_relevant)
            returnParameter!(q,X_mat,
                            Z,X_0_t,X_t,X_last,X_int,draws,F_t,
                            γlen,Ilen,β0len,βlen,σlen)
            @inbounds Y_list[q_i] = X_mat[:,:]

            if (γlen < q <= Ilen)
                grad_calc_ω!(dω_x,ω_i,X_mat)


                # dω_list[q_i] = dω_x[:]
                # dS_u_x[:] .= 0.0
                # dS_uncond_list[q_i] = dS_u_x[:]

                # @inbounds @fastmath @simd for k in 1:K
                #     dS_x[k] = dω_x[k]*(s_uncond[k] - y_last[k])
                # end
                dS_x[:] = dω_x.*s_adjust
            else
                grad_calc_s!(dS_u_x,s_n,
                            μ_ij,X_mat,
                            μ_ij_sums)

                # dS_uncond_list[q_i] = dS_u_x[:]
                # dω_x[:] .= 0.0
                # dω_list[q_i] = dω_x[:]

                dS_x[:] = (dS_u_x[:].*ω_i)./N
            end

            ## Calculate Gradient
            grad[q]+= combine_grad(gll_t1,dS_x)
        end

    return ll_obs,pars_relevant
end


function combine_hess(gll_t1::Vector{T},gll_t2::Vector{T},
                    dS_xy::Vector{T},dS_x::Vector{T},
                    dS_y::Vector{T}) where T
    hess_obs = 0.0
    K = length(dS_xy)

    @inbounds @fastmath @simd for k in 1:K
        hess_obs += gll_t1[k]*(dS_xy[k]) - gll_t2[k]*(dS_x[k])*(dS_y[k])
    end
    return hess_obs
end


function combine_grad(gll_t1::Vector{T},
                    dS_x::Vector{T}) where T
    grad_obs = 0.0
    K = length(dS_x)

    @inbounds @fastmath @simd for k in 1:K
        grad_obs += gll_t1[k]*dS_x[k]
    end
    return grad_obs
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    X_mat::Matrix{Float64},Y::Matrix{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]
        @inbounds @fastmath y = Y[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_x!(n::Int64,s_n::Vector{T},
                    X_mat::Matrix{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
    end

    return dμ_ij_x_sums
end


function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Matrix{Float64},Y::Matrix{Float64},
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat[n,k] - Γ_x)*(Y[n,k] - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end


function calc_prodTerms_x!(n::Int64,
                        dS_x::Vector{T},
                        X_mat::Matrix{Float64},
                        s_n::Vector{T},
                        Γ_x::T) where T

    K = length(dS_x)
    for k in 1:K
        @inbounds @fastmath tx = s_n[k]*(X_mat[n,k] - Γ_x)
        @inbounds @fastmath dS_x[k]+= tx
    end
end


function hess_calc!(dS_xy::Vector{Float64},s_n::Vector{Float64},
                    μ_ij::Matrix{Float64},
                    X_mat::Matrix{Float64},Y::Matrix{Float64},
                    μ_ij_sums::Vector{Float64})
    dS_xy[:] .= 0.0


    (N,K) = size(μ_ij)

    ### 0 Draw Calculation ###
    for n in 1:N
        μ_ij_sums_n = μ_ij_sums[n]

        dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums = calc_derSums_xy!(n,s_n,X_mat,Y,μ_ij,μ_ij_sums_n)

        @fastmath Γ_x = dμ_ij_x_sums/μ_ij_sums_n
        @fastmath Γ_y = dμ_ij_y_sums/μ_ij_sums_n
        @fastmath Γ_xy = dμ_ij_xy_sums/μ_ij_sums_n
        # s_n = μ_ij[n,:].*δ/μ_ij_sums_n

        calc_prodTerms_xy!(n,dS_xy,
                            X_mat,Y,s_n,
                            Γ_x,Γ_y,Γ_xy)
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Matrix{Float64},Y_mat::Matrix{Float64},
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat[1,k]*Y_mat[1,k]*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function grad_calc_s!(dS_x::Vector{T},
                    s_n::Vector{T},
                    μ_ij::Matrix{T},
                    X_mat::Matrix{Float64},
                    μ_ij_sums::Vector{T}) where T

    dS_x[:] .= 0.0

    (N,K) = size(μ_ij)

    for n in 1:N
        μ_ij_sums_n = μ_ij_sums[n]
        dμ_ij_x_sums = calc_derSums_x!(n,s_n,X_mat,μ_ij,μ_ij_sums_n)

        @fastmath Γ_x = dμ_ij_x_sums/μ_ij_sums_n
        # s_n = s_n./μ_ij_sums_n

        calc_prodTerms_x!(n,dS_x,X_mat,s_n,Γ_x)
    end
    return nothing
end

function grad_calc_ω!(dω_x::Vector{T},
                    ω::T,
                    X_mat::Matrix{Float64}) where T

    dω_x[:] .= 0.0
    K = length(dω_x)
    for k in 1:K
        dω_x[k] = (X_mat[1,k]*ω*(1-ω))#*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function getIndex_broad(A::Float64,n::Int64,k::Int64)
    return A
end

function getIndex_broad(A::Array{Float64,1},n::Int64,k::Int64)
    return A[k]
end

function getIndex_broad(A::Array{Float64,2},n::Int64,k::Int64)
    return A[n,k]
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

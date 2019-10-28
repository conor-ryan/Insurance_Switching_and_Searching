using NLopt
using ForwardDiff


function estimate_ng!(d::InsuranceLogit, p0;method=:LN_SBPLX)
    # Set up the optimization
    opt = Opt(method, length(p0))
    ftol_rel!(opt,1e-8)
    xtol_rel!(opt,1e-8)
    maxtime!(opt, 500000)

    lb = repeat([-100],inner=length(p0))
    ub = repeat([100],inner=length(p0))

    # lower_bounds!(opt_stage1, lb)
    # upper_bounds!(opt_stage1, ub)

    ll(x) = log_likelihood(d,x)
    disp_length = min(20,length(p0))
    count = 0
    function ll(x, grad)
        count +=1
        x_displ = x[1:disp_length]
        println("Iteration $count at $x_displ")
        obj = ll(x)
        println("Objective equals $obj on iteration $count")

        return obj
    end

    # Set Objective
    max_objective!(opt, ll)

    # Run Optimization
    minf, minx, ret = optimize(opt, p0)
    println("Got $minf at $minx after $count iterations (returned $ret)")

    # Return the object
    return ret, minf, minx
end


function newton_raphson_ll(d,p0;grad_tol=1e-8,f_tol=1e-8,x_tol=1e-8,
    max_itr=2000,strict=true,Hess_Skip_Steps=10)
    ## Initialize Parameter Vector
    p_vec = p0
    N = length(p0)

    cnt = 0
    grad_size = 10000
    f_eval_old = 1.0
    # # Initialize δ
    param_dict = parDict(d,p_vec)
    Pop = length(d.data._personIDs)

    # Initialize Gradient
    grad_new = similar(p0)
    hess_new = Matrix{Float64}(undef,length(p0),length(p0))
    Eye = Matrix{Float64}(1.0I,length(p0),length(p0))
    f_final_val = 0.0
    max_trial_cnt = 0
    NaN_steps = 0
    trial_end = 5
    hess_steps = 0
    p_last = copy(p_vec)
    grad_last = copy(grad_new)
    H_last = copy(hess_new)
    disp_length = min(length(p0),20)
    f_min = -1e3
    p_min  = similar(p_vec)
    no_progress=0
    flag = "empty"
    if strict
        mistake_thresh = 1.00
    else
        mistake_thresh = 1.25
    end

    ## Initialize Step
    step = 1
    real_hessian=0

    ## Bound Index
    bound_ind =  Int.((length(p0) - d.parLength[:FE] - d.parLength[:σ] + 1):((length(p0) - d.parLength[:FE])))
    if d.parLength[:σ]==0
        bound_ind = Vector{Int64}(undef,0)
    end
    println("Bounded Parameters: $bound_ind")
    constrained = 0
    constraint = -1e4

    ## Tolerance Counts
    f_tol_cnt = 0
    x_tol_cnt = 0
    ga_conv_cnt = 0
    skip_x_tol = 0
    ga_cnt = 0
    ga_max_itr = 10
    # Maximize by Newtons Method
    while (grad_size>grad_tol) & (cnt<max_itr) & (max_trial_cnt<20)
        cnt+=1
        trial_cnt=0

        ## Check Constraint
        if any(p_vec[bound_ind].<constraint)
            ind = bound_ind[findall(p_vec[bound_ind].<constraint)]
            println("Hit Constraint at $ind")
            p_vec[ind] = abs.(p_vec[ind])
            if any(abs.(p_vec[bound_ind]).<constraint)
                ind = bound_ind[findall(p_vec[bound_ind].<constraint)]
                println("Hit Constraint at $ind, in magnitude")
                p_vec[ind].= constraint
                constrained = 1
            else
                constrained = 0
            end
        else
            constrained = 0
        end

        # Compute Gradient, holding δ fixed
        if hess_steps==0
            println("Compute Hessian")
            fval = log_likelihood!(hess_new,grad_new,d,p_vec)
            update, H_k = boundAtZero(bound_ind,p_vec,hess_new,grad_new,constraint)
            # H_k = inv(hess_new)
            real_hessian=1
            if cnt==1
                f_min = fval
                p_min = p_vec
                println("Initial Minimum Value: $f_min")
            end
        else
            println("BFGS Approximation")
            fval = log_likelihood!(grad_new,d,p_vec)
            Δxk = p_vec - p_last
            yk = grad_new - grad_last
            # H_k = (Eye - (Δxk*yk')./(yk'*Δxk) )*H_last*(Eye - (yk*Δxk')./(yk'*Δxk) ) + (Δxk*Δxk')./(yk'*Δxk)
            update, H_k = boundAtZero(bound_ind,p_vec,Eye,H_last,Δxk,yk,grad_new,constraint)
            real_hessian=0
        end
        # H_k = inv(hess_new)
        # update = -H_k*grad_new

        if (fval>f_min) | (constrained == 1)
            if (abs(fval-f_min)<f_tol) & (skip_x_tol==0)
                f_tol_cnt += 1
            end
            if (maximum(abs.(p_vec - p_min))<x_tol) & (skip_x_tol==0)
                x_tol_cnt += 1
            end

            f_min = copy(fval)
            p_min[:] = p_vec[:]

            no_progress=0
        else
            no_progress+=1
        end
        skip_x_tol = 0


        grad_size = maximum(abs.(grad_new))
        if (grad_size<grad_tol) |(f_tol_cnt>1) | (x_tol_cnt>1) | (ga_conv_cnt>1)
            println("Got to Break Point...?")
            println(grad_size)
            println(f_tol_cnt)
            println(x_tol_cnt)
            println(ga_conv_cnt)
            flag = "converged"
            break
        end
        if strict==false
            if grad_size>.1
                mistake_thresh = 1.25
            else
                mistake_thresh = 1.05
            end
        end


        if any(isnan.(update))
            println("Step contains NaN")
            println("Algorithm Failed")
            return p_min,f_min

            # #Check Hessian
            # eig = sort(abs.(eigvals(hess_new)))
            # sm_e = eig[1]
            # println("Smallest Eigenvalue: $sm_e ")
            # NaN_steps +=1
            # grad_size = sqrt(dot(grad_new,grad_new))
            # update = -(1/grad_size).*grad_new
        else
            NaN_steps = 0
        end


        if no_progress>5
            no_progress = 0
            println("Return: Limit on No Progress")
            p_vec = copy(p_min)
            fval = log_likelihood!(grad_new,d,p_vec)
            grad_size = maximum(abs.(grad_new))
            step = 1/grad_size
            update = - step.*grad_new
            mistake_thresh = 1.00
        end

        step_size = maximum(abs.(update))
        if step_size>10
        update = update./step_size
        ind = findall(abs.(update).==1)
        val_disp = p_vec[ind]
            println("Max Parameter Adjustment: $ind, $val_disp")
        step_size = 1
        end


        p_test = p_vec .+ update
        f_test = log_likelihood(d,p_test)

        if hess_steps<Hess_Skip_Steps
            hess_steps+=1
        else
            hess_steps=0
        end

        if ga_cnt>=2
            ga_max_itr = 50
        else
            ga_max_itr = 10
        end

        trial_max = 0
        while ((f_test<fval*mistake_thresh) | isnan(f_test)) & (trial_max==0)
            if real_hessian==0
                skip_x_tol = 1
            end
            if trial_cnt==0
                p_test_disp = p_test[vcat(1:10,bound_ind)]
                println("Trial (Init): Got $f_test at parameters $p_test_disp")
                println("Previous Iteration at $fval")
            end
            if trial_cnt<=2
                update/= 20
            else
                update/= 200
            end
            step_size = maximum(abs.(update))
            if (step_size>x_tol)
                p_test = p_vec .+ update
                f_test = log_likelihood(d,p_test)
                p_test_disp = p_test[vcat(1:10,bound_ind)]
                println("Trial (NR): Got $f_test at parameters $p_test_disp")
                println("Previous Iteration at $fval")
                trial_cnt+=1
            elseif real_hessian==1
                hess_steps = 0
                trial_max = 1
                ga_cnt+=1
                println("RUN ROUND OF GRADIENT ASCENT")
                p_test, f_test, ga_conv = gradient_ascent(d,p_vec,max_itr=ga_max_itr,strict=true)
                ga_conv_cnt += ga_conv
            else
                println("No Advancement")
                hess_steps = 0
                p_test = copy(p_vec)
                break
            end
        end
        if (step_size>x_tol)
            ga_cnt = 0
        end
        if NaN_steps>5
            println("Hessian might be singular")
            println("RUN ROUND OF GRADIENT ASCENT")
            p_test, f_test,ga_conv = gradient_ascent(d,p_test,max_itr=20,strict=true)
            ga_conv_cnt += ga_conv
        end

        p_last = copy(p_vec)
        p_vec = copy(p_test)
        grad_last = copy(grad_new)
        H_last = copy(H_k)
        p_vec_disp = p_vec[vcat(1:10,bound_ind)]
        f_final_val = fval
        println("Update Parameters to $p_vec_disp")


        println("Gradient Size: $grad_size")
        println("Step Size: $step_size")
        f_full = f_min*Pop
        println("Function Value is $f_test, $f_full at iteration $cnt")
        println("Steps since last improvement: $no_progress")
    end
    f_full = f_min*Pop
    println("Lowest Function Value is $f_min, $f_full at $p_min")
    return p_min,f_min
end

#
# function newton_raphson_ll(d,p0;grad_tol=1e-8,f_tol=1e-8,x_tol=1e-8,step_tol=1e-8,max_itr=2000)
#     ## Initialize Parameter Vector
#     p_vec = p0
#     N = length(p0)
#
#     Pop = sum(weight(d.data).*choice(d.data))
#
#     cnt = 0
#     grad_size = 10000
#     f_eval_old = 1.0
#     # # Initialize δ
#     param_dict = parDict(d,p_vec)
#
#     disp_length = min(length(p0),20)
#     # Initialize Gradient
#     grad_new = similar(p0)
#     hess_new = Matrix{Float64}(undef,length(p0),length(p0))
#     f_final_val = 0.0
#     max_trial_cnt = 0
#     flag = "empty"
#     f_max = -1e3
#     p_min  = similar(p_vec)
#     no_progress = 0
#     f_tol_cnt = 0
#     x_tol_cnt = 0
#     mistake_thresh=1.00
#     # Maximize by Newtons Method
#     while (grad_size>grad_tol) & (cnt<max_itr) & (max_trial_cnt<20)
#         cnt+=1
#         trial_cnt = 0
#
#         # Compute Gradient, holding δ fixed
#
#         fval = log_likelihood!(hess_new,grad_new,d,p_vec)
#         if cnt==1
#             fval_pop = fval*Pop
#             println("Function Value is $fval_pop at iteration 0")
#         end
#         if (cnt==1) | (fval>f_max)
#             if abs(fval-f_max)<f_tol
#                 f_tol_cnt += 1
#             end
#             if maximum(abs.((p_vec - p_min)./p_min))<x_tol
#                 x_tol_cnt += 1
#             end
#
#             f_max = copy(fval)
#             p_min[:] = p_vec[:]
#
#             no_progress=0
#         else
#             f_tol_cnt = 0
#             x_tol_cnt = 0
#             no_progress+=1
#         end
#
#         grad_size = sqrt(dot(grad_new,grad_new))
#         if (grad_size<grad_tol) |(f_tol_cnt>3) | (x_tol_cnt>3)
#             println("Got to Break Point")
#             println(grad_size)
#             println(f_tol_cnt)
#             println(x_tol_cnt)
#             flag = "converged"
#             break
#         end
#
#
#
#         # ForwardDiff.gradient!(grad_new, ll, p_vec)
#         # println("Gradient is $grad_new")
#         #
#         #
#         # hess_new = Matrix{Float64}(N,N)
#         # ForwardDiff.hessian!(hess_new, ll, p_vec)
#         # println("Hessian is $hess_new")
#         # if any(abs.(diag(hess_new)).<1e-10)
#         #     return p_vec, hess_new
#         # end
#
#         update = -inv(hess_new)*grad_new
#         if any(isnan.(update))
#             println("Step contains NaN")
#             #Check Hessian
#             # eig = sort(abs.(eigvals(hess_new)))
#             # sm_e = eig[1]
#             # println("Smallest Eigenvalue: $sm_e ")
#             NaN_steps +=1
#             grad_size = sqrt(dot(grad_new,grad_new))
#             update = -(1/grad_size).*grad_new
#         else
#             NaN_steps = 0
#         end
#
#         step_size = maximum(abs.(update))
#         if step_size>20
#             ind = findall(abs.(update).>20)
#             if any(abs.(p_vec[ind]).>200)
#                 update = (update./step_size).*5
#                 step_size = 5
#             else
#                 update = (update./step_size).*20
#                 step_size = 20
#             end
#             ind = findall(abs.(update).==step_size)
#             val_disp = p_vec[ind]
#             println("Max Parameter Adjustment: $ind, $val_disp")
#         end
#
#         if no_progress>5
#             flag = "no better point"
#             println("No Progress in Algorithm")
#             p_test, f_test = gradient_ascent(d,p_vec,max_itr=50,strict=false)
#         else
#             p_test = p_vec .+ update
#             f_test = log_likelihood(d,p_test)
#         end
#
#         step_size_thresh = minimum(vcat(step_size,maximum(abs.(update./p_vec))))
#         trial_max = 0
#         while ((f_test<fval*mistake_thresh) | isnan(f_test)) & (trial_max==0)
#             if trial_cnt==0
#                 p_test_disp = p_test[1:20]
#                 println("Trial (Init): Got $f_test at parameters $p_test_disp")
#                 println("Previous Iteration at $fval")
#             end
#             if (step_size_thresh>x_tol)
#                 if trial_cnt<=2
#                     update/= 20
#                 else
#                     update/= 200
#                 end
#                 step_size = maximum(abs.(update))
#                 step_size_thresh = minimum(vcat(step_size,maximum(abs.(update./p_vec))))
#                 p_test = p_vec .+ update
#                 f_test = log_likelihood(d,p_test)
#                 p_test_disp = p_test[1:20]
#                 println("Trial (NR): Got $f_test at parameters $p_test_disp")
#                 println("Previous Iteration at $fval")
#                 trial_cnt+=1
#             else
#                 trial_max = 1
#                 println("RUN ROUND OF GRADIENT ASCENT")
#                 p_test, f_test = gradient_ascent(d,p_vec,max_itr=5,strict=true)
#             end
#         end
#         p_vec = copy(p_test)
#         p_vec_disp = p_vec[1:disp_length]
#         f_final_val = f_test
#         println("Update Parameters to $p_vec_disp")
#
#         fval_pop = f_test*Pop
#
#
#         println("Gradient Size: $grad_size")
#         println("Function Value is $f_test, $fval_pop at iteration $cnt")
#     end
#     # if (grad_size>grad_tol)
#     #     println("Estimate Instead")
#     #     ret, f_final_val, p_vec = estimate!(d,p0)
#
#
#     return p_min,f_max,flag
# end



function gradient_ascent(d,p0;grad_tol=1e-8,f_tol=1e-8,x_tol=1e-8,max_itr=2000,strict=false)
    ## Initialize Parameter Vector
    p_vec = copy(p0)
    N = length(p0)

    cnt = 0
    grad_size = 10000
    f_eval_old = 1.0
    # # Initialize δ
    param_dict = parDict(d,p_vec)

    # Initialize Gradient
    grad_new = similar(p0)
    hess_new = Matrix{Float64}(undef,length(p0),length(p0))
    f_final_val = 0.0
    max_trial_cnt = 0
    p_last = copy(p_vec)
    grad_last = copy(grad_new)
    disp_length = min(length(p0),20)
    f_max = -1e3
    p_min  = similar(p_vec)
    no_progress=0
    mistake_thresh = 1.25
    restart = 0
    flag = 0

    ## Bound Index
    bound_ind =  Int.((length(p0) - d.parLength[:FE] - d.parLength[:σ] + 1):((length(p0) - d.parLength[:FE])))
    if d.parLength[:σ]==0
        bound_ind = Vector{Int64}(undef,0)
    end
    println("Bounded Parameters: $bound_ind")
    constrained = 0
    constraint = -1e4

    init_val = log_likelihood(d,p_vec)
    println("Initial Value is $init_val")

    ### Tolerance Counts
    f_tol_cnt = 0
    x_tol_cnt = 0

    real_gradient = 1
    grad_steps = 0
    Grad_Skip_Steps = 15
    step = 1.0
    fval = -1e3
    # Maximize by Newtons Method
    while (grad_size>grad_tol) & (cnt<max_itr) & (max_trial_cnt<20)
        cnt+=1
        trial_cnt = 0

        ## Check Constraint
        if any(p_vec[bound_ind].<constraint)
            ind = bound_ind[findall(p_vec[bound_ind].<constraint)]
            println("Hit Constraint at $ind")
            p_vec[ind] = abs.(p_vec[ind])
            if any(abs.(p_vec[bound_ind]).<constraint)
                ind = bound_ind[findall(p_vec[bound_ind].<constraint)]
                println("Hit Constraint at $ind, in magnitude")
                p_vec[ind].= constraint
                constrained = 1
            else
                constrained = 0
            end
        else
            constrained = 0
        end

        # Compute Gradient, holding δ fixed


        if grad_steps==0
            println("Compute Gradient")
            fval = log_likelihood!(grad_new,d,p_vec)
            real_gradient=1
            if cnt==1
                f_max = fval
                p_min = p_vec
                println("Initial Minimum Value: $f_max")
            end
        else
            # fval = GMM_objective(d,p_vec,W)
            real_gradient=0
        end

        # grad_size = sqrt(dot(grad_new,grad_new))
        grad_size = maximum(abs.(grad_new))
        if grad_size<.1
            mistake_thresh = 1.05
        else
            mistake_thresh = 1.25
        end
        if strict
            mistake_thresh = 1.00
        end

        if (grad_size<grad_tol) | (x_tol_cnt>1) |( f_tol_cnt>1)
            println("Got to Break Point...?")
            println(grad_size)
            println(x_tol_cnt)
            println(f_tol_cnt)
            flag = 1
            break
        end

        if cnt==1
            step = min(1/grad_size,1.0)
        elseif (real_gradient==1)
            g = p_vec - p_last
            y = grad_new - grad_last
            step = abs(dot(g,g)/dot(g,y))
        end

        if restart == 1
            flag = "no progress"
            step = min(1/grad_size,1.0)
            restart = 0
            no_progress = 0
            mistake_thresh = 1.00
        end

        update = boundAtZero(bound_ind,p_vec,step,grad_new,constraint)

        p_test = p_vec .+ update

        f_test = log_likelihood(d,p_test)

        if (f_test<fval) & (real_gradient==0)
            grad_steps=0
            continue
        elseif (grad_steps<Grad_Skip_Steps)
            grad_steps+=1
        else
            grad_steps=0
        end

        while ((f_test<fval*mistake_thresh) | isnan(f_test)) & (step>1e-15)
            p_test_disp = p_test[1:disp_length]
            if trial_cnt==0
                println("Trial: Got $f_test at parameters $p_test_disp")
                println("Previous Iteration at $fval")
                println("Reducing Step Size...")
            end
            update/= 20
            p_test = p_vec .+ update
            f_test = log_likelihood(d,p_test)
                println("Trial: Got $f_test at parameters $p_test_disp")
                println("Previous Iteration at $fval")
            trial_cnt+=1
            if (maximum(abs.(update))<x_tol) & (real_gradient==0)
                println("Failed with Approximate Gradient")
                break
            elseif (maximum(abs.(update))<1e-16)
                println("Failed Step is too small")
                break
            end
        end

        ## Update Minimum Value
        if (f_test>f_max) | (constrained==1)
            if (abs(f_test-f_max)<f_tol) & (real_gradient==1)
                f_tol_cnt += 1
            end
            if (maximum(abs.(p_test - p_min))<x_tol) & (real_gradient==1)
                    x_tol_cnt += 1
            end
            f_max = copy(f_test)
            p_min[:] = p_test[:]
            no_progress=0
        else
            no_progress+=1
            grad_steps=0
        end

        if real_gradient==1
            p_last = copy(p_vec)
            grad_last = copy(grad_new)
        end
        p_vec = copy(p_test)
        fval = copy(f_test)
        p_vec_disp = p_vec[vcat(1:10,bound_ind)]
        f_final_val = fval
        println("Update Parameters to $p_vec_disp")
        println("Gradient Size: $grad_size")
        println("Step Size: $step")
        println("Function Value is $fval at iteration $cnt")
        println("Steps since last improvement: $no_progress")


        if (no_progress==15)
            println("No Progress - Restart near p_min")
            println("Lowest Function Value so far: $f_max")
            restart = 1
            p_vec = copy(p_min)
            # flag = "no progress"
            # break
        end
    end

    println("Lowest Function Value is $f_max at $p_min")
    return p_min,f_max, flag
end

function boundAtZero(p_ind::Vector{Int64},p::Vector{Float64},
    α::Float64,grad::Vector{Float64},bound::Float64)

    update = α*grad
    if length(p_ind)==0
        return update
    end

    bound_ind = p_ind[findall((p[p_ind].<=bound) .& (update[p_ind].<0))]

    if length(bound_ind)==0
        return update
    end

    update = zeros(length(p))
    unbounded = .!(inlist(1:length(p),bound_ind))

    update[unbounded] = α*grad[unbounded]
    return update
end

function boundAtZero(p_ind::Vector{Int64},p::Vector{Float64},hess::Matrix{Float64},grad::Vector{Float64},bound::Float64)
    check = issuccess(cholesky(-hess,check=false))
    if !check
        println("Not Concave: Compute Direction of Sufficient Ascent") #http://web.stanford.edu/class/cme304/docs/newton-type-methods.pdf
        E = eigen(hess)
        max_eig_val = maximum(E.values)
        println("Max Eigenvalue is $max_eig_val")
        Λ = -abs.(Diagonal(E.values))
        hess = E.vectors*Λ*E.vectors'

        if minimum(abs.(E.values))<1e-10
            println("An Eigenvalue is close to 0")
            max_param = maximum(abs.(p))
            max_index = findall(abs.(p).==max_param)
            println("Parameter w/ largest magnitude: $max_index is $max_param")
        end
        # (J,K) = size(H_adj)
        # for j in 1:J
        #     for k in 1:j
        #         H_adj[j,k] = H_adj[k,j]
        #     end
        # end
    end
    H_k_f = inv(hess)
    update = -H_k_f*grad

    if length(p_ind)==0
        println("Standard Update 1")
        return update, H_k_f
    end

    bound_ind = p_ind[findall((p[p_ind].<=bound) .& (update[p_ind].<0))]

    if length(bound_ind)==0
        println("Standard Update 2")
        return update, H_k_f
    end

    update = zeros(length(p))
    unbounded = .!(inlist(1:length(p),bound_ind))

    H_k = inv(hess[unbounded,unbounded])
    update[unbounded] = -H_k*grad[unbounded]

    H_k_f[unbounded,unbounded] = H_k

    println("Bounded Update")
    return update, H_k_f
end

function boundAtZero(p_ind::Vector{Int64},p::Vector{Float64},
    Eye::Matrix{Float64},H_last::Matrix{Float64},Δx::Vector{Float64},Δy::Vector{Float64},
    grad::Vector{Float64},bound::Float64)
    H_k_f = (Eye - (Δx*Δy')./(Δy'*Δx) )*H_last*(Eye - (Δy*Δx')./(Δy'*Δx) ) + (Δx*Δx')./(Δy'*Δx)
    update = -H_k_f*grad

    if length(p_ind)==0
        println("Standard Update 1")
        return update, H_k_f
    end

    bound_ind = p_ind[findall((p[p_ind].<=bound) .& (update[p_ind].<0))]

    if length(bound_ind)==0
        println("Standard Update 2")
        return update, H_k_f
    end

    update = zeros(length(p))
    unbounded = .!(inlist(1:length(p),bound_ind))

    Δxk = Δx[unbounded]
    Δyk = Δy[unbounded]

    H_k = (Eye[unbounded,unbounded] - (Δxk*Δyk')./(Δyk'*Δxk) )*H_last[unbounded,unbounded]*(Eye[unbounded,unbounded] - (Δyk*Δxk')./(Δyk'*Δxk) ) + (Δxk*Δxk')./(Δyk'*Δxk)
    update[unbounded] = -H_k*grad[unbounded]

    H_k_f[unbounded,unbounded] = H_k

    println("Bounded Update")
    return update, H_k_f
end

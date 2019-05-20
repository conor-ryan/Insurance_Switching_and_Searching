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


function newton_raphson_ll(d,p0;grad_tol=1e-8,f_tol=1e-8,x_tol=1e-8,step_tol=1e-8,max_itr=2000)
    ## Initialize Parameter Vector
    p_vec = p0
    N = length(p0)

    Pop = sum(weight(d.data).*choice(d.data))

    cnt = 0
    grad_size = 10000
    f_eval_old = 1.0
    # # Initialize δ
    param_dict = parDict(d,p_vec)

    disp_length = min(length(p0),20)
    # Initialize Gradient
    grad_new = similar(p0)
    hess_new = Matrix{Float64}(undef,length(p0),length(p0))
    f_final_val = 0.0
    max_trial_cnt = 0
    flag = "empty"
    f_max = -1e3
    p_min  = similar(p_vec)
    no_progress = 0
    f_tol_cnt = 0
    x_tol_cnt = 0
    mistake_thresh=1.00
    # Maximize by Newtons Method
    while (grad_size>grad_tol) & (cnt<max_itr) & (max_trial_cnt<20)
        cnt+=1
        trial_cnt = 0

        # Compute Gradient, holding δ fixed

        fval = log_likelihood!(hess_new,grad_new,d,p_vec)
        if cnt==1
            fval_pop = fval*Pop
            println("Function Value is $fval_pop at iteration 0")
        end
        if (cnt==1) | (fval>f_max)
            if abs(fval-f_max)<f_tol
                f_tol_cnt += 1
            end
            if maximum(abs.((p_vec - p_min)./p_min))<x_tol
                x_tol_cnt += 1
            end

            f_max = copy(fval)
            p_min[:] = p_vec[:]

            no_progress=0
        else
            f_tol_cnt = 0
            x_tol_cnt = 0
            no_progress+=1
        end

        grad_size = sqrt(dot(grad_new,grad_new))
        if (grad_size<grad_tol) |(f_tol_cnt>3) | (x_tol_cnt>3)
            println("Got to Break Point")
            println(grad_size)
            println(f_tol_cnt)
            println(x_tol_cnt)
            flag = "converged"
            break
        end



        # ForwardDiff.gradient!(grad_new, ll, p_vec)
        # println("Gradient is $grad_new")
        #
        #
        # hess_new = Matrix{Float64}(N,N)
        # ForwardDiff.hessian!(hess_new, ll, p_vec)
        # println("Hessian is $hess_new")
        # if any(abs.(diag(hess_new)).<1e-10)
        #     return p_vec, hess_new
        # end

        update = -inv(hess_new)*grad_new
        if any(isnan.(update))
            println("Step contains NaN")
            #Check Hessian
            eig = sort(abs.(eigvals(hess_new)))
            sm_e = eig[1]
            println("Smallest Eigenvalue: $sm_e ")
            NaN_steps +=1
            grad_size = sqrt(dot(grad_new,grad_new))
            update = -(1/grad_size).*grad_new
        else
            NaN_steps = 0
        end

        if no_progress>5
            flag = "no better point"
            println("No Progress in Algorithm")
            p_test, f_test = gradient_ascent(d,p_vec,max_itr=50,strict=false)
        else
            p_test = p_vec .+ update
            f_test = log_likelihood(d,p_test)
        end





        step_size = maximum(abs.(update))
        step_size_thresh = minimum(vcat(step_size,maximum(abs.(update./p_vec))))
        trial_max = 0
        while ((f_test<fval*mistake_thresh) | isnan(f_test)) & (trial_max==0)
            if trial_cnt==0
                p_test_disp = p_test[1:20]
                println("Trial (Init): Got $f_test at parameters $p_test_disp")
                println("Previous Iteration at $fval")
            end
            if (step_size_thresh>x_tol)
                if trial_cnt<=2
                    update/= 20
                else
                    update/= 200
                end
                step_size = maximum(abs.(update))
                step_size_thresh = minimum(vcat(step_size,maximum(abs.(update./p_vec))))
                p_test = p_vec .+ update
                f_test = log_likelihood(d,p_test)
                p_test_disp = p_test[1:20]
                println("Trial (NR): Got $f_test at parameters $p_test_disp")
                println("Previous Iteration at $fval")
                trial_cnt+=1
            else
                trial_max = 1
                println("RUN ROUND OF GRADIENT ASCENT")
                p_test, f_test = gradient_ascent(d,p_vec,max_itr=5,strict=true)
            end
        end
        p_vec = copy(p_test)
        p_vec_disp = p_vec[1:disp_length]
        f_final_val = f_test
        println("Update Parameters to $p_vec_disp")

        fval_pop = f_test*Pop


        println("Gradient Size: $grad_size")
        println("Function Value is $f_test, $fval_pop at iteration $cnt")
    end
    # if (grad_size>grad_tol)
    #     println("Estimate Instead")
    #     ret, f_final_val, p_vec = estimate!(d,p0)


    return p_min,f_max,flag
end



function gradient_ascent(d,p0;grad_tol=1e-8,f_tol=1e-8,x_tol=1e-8,max_itr=2000,strict=false)
    ## Initialize Parameter Vector
    p_vec = p0
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
    flag = "empty"

    ### Tolerance Counts
    f_tol_cnt = 0
    x_tol_cnt = 0

    # Maximize by Newtons Method
    while (grad_size>grad_tol) & (cnt<max_itr) & (max_trial_cnt<20)
        cnt+=1
        trial_cnt = 0

        # Compute Gradient, holding δ fixed

        fval = log_likelihood!(grad_new,d,p_vec)
        if (cnt==1) | (fval>f_max)
            if abs(fval-f_max)<f_tol
                f_tol_cnt += 1
            end
            if maximum(abs.(p_vec - p_min))<x_tol
                x_tol_cnt += 1
            end

            f_max = copy(fval)
            p_min[:] = p_vec[:]

            no_progress=0
        else
            no_progress+=1
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
            flag = "near zero"
            break
        end

        if cnt==1
            step = 1/grad_size
        else
            g = p_vec - p_last
            y = grad_new - grad_last
            step = abs(dot(g,g)/dot(g,y))
        end

        if restart == 1
            flag = "no progress"
            step = (1/grad_size)
            restart = 0
            no_progress = 0
            mistake_thresh = 1.00
        end

        p_test = p_vec .+ step.*grad_new

        f_test = log_likelihood(d,p_test)

        while ((f_test<fval*mistake_thresh) | isnan(f_test)) & (trial_cnt<10)
            p_test_disp = p_test[1:disp_length]
            if trial_cnt==0
                println("Trial: Got $f_test at parameters $p_test_disp")
                println("Previous Iteration at $fval")
                println("Reducing Step Size...")
            end
            step/= 2
            p_test = p_vec .+ step.*grad_new
            f_test = log_likelihood(d,p_test)
            trial_cnt+=1
            if (trial_cnt==10) & (grad_size>1e-5)
                println("Algorithm Stalled: Random Step")
                max_trial_cnt+=1
                step = rand(length(step))/1000 .-.005
            elseif (trial_cnt==10) & (grad_size<=1e-5)
                println("Algorithm Stalled: Random Step")
                max_trial_cnt+=1
                step = rand(length(step))/10000 .-.005
            end
        end

        p_last = copy(p_vec)
        p_vec = copy(p_test)
        grad_last = copy(grad_new)
        p_vec_disp = p_vec[1:20]
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

using JLD
using PyPlot


# % constant 
    WC =     [2.285714286 0.1428571429 -0.4285714286 ; 
              0.1428571429 1.571428571 0.2857142857;
              -0.4285714286 0.2857142857 1.142857143];

    # % linear 
    WL =   [0 -4.57142857142857 0 ; 
              -4.57142857142857 -0.571549759244019 0.857142857142854;
              0 0.857142857142854 0];

    # % quadratic 
    WQ =    [0 0 0 ; 
              0 9.14297833067258 0;
              0 0 0]; # setting up interpolating nodes CGL



ps_params(config) = Dict("nodes" => CGLnodes(config),
                         "weights" => ClenshawCurtisWeight(config)
                         )
sparsify(x, eps) = abs(x) < eps ? 0.0 : x
@vectorize_2arg Float64 sparsify

CGLnodes(config) = [0.5 - 0.5*cos(k*pi/config["N"]) for k in 0:config["N"]]

LagrangePolynomial(config, ps, t) = 
sparsify(prod(hcat([ [(t - ps["nodes"][i])/(ps["nodes"][j] - ps["nodes"][i])   for i in 1:config["N"]+1 if i!=j] for j in 1:config["N"]+1]...),1), config["sparse_eps"])'

# outputs: vector D(t_i) where each element is D_j(t_i) - the deriative of the jth polynomial evaluated at point t_i
DiffLagrangePolynomial(config, ps, t) = sparsify([sum([prod([(t - ps["nodes"][m])/(ps["nodes"][j] - ps["nodes"][m]) for m in 1:config["N"]+1 if m != j && m != n])/(ps["nodes"][j] - ps["nodes"][n]) for n in 1:config["N"]+1 if n!=j]) for j in 1:config["N"]+1], config["sparse_eps"])

function  ChebyshevPolynomial(config, s)
    s = 2*s - 1
    T0 = ones(length(s))
    T1 = s
    T = hcat(T0, T1)
    if config["deg"] >= 2
        for i = 2:config["deg"]
            T2 = 2*s.*T1 - T0;
            T =  hcat(T, T2)
            T0 = T1;
            T1 = T2;
        end
    end
    T'
    # nodes go across
    # basis goes down
end

function ChebyshevSecondKind(config, s, k = 0)
    deg = config["deg"] - k
    s = 2*s - 1
    U0 = ones(length(s))
    U1 =  2*s
    if deg == 0
        U = U0;
    elseif deg == 1
        U = hcat(U0, U1)
    else 
        U = hcat(U0, U1)
    end
    if deg >= 2
        for i = 2:deg
            U2 = 2*s*U1 - U0
            U =  hcat(U, U2)
            U0 = U1
            U1 = U2
        end
    end
    U'
end

function DiffChebyshevPolynomial(config, s)
    bot = vcat(0, ChebyshevSecondKind(config, s, 1))
    sparsify(2*collect(0:config["deg"]).*bot, config["sparse_eps"])
end


function Constraints(config, ps, xstar, xcurr)
    Z = zeros(config["deg"] + 1)'
    start = 0
    finish = 0
    if config["type"] == "lagrange"
        start = LagrangePolynomial(config, ps, 0)'
        finish = LagrangePolynomial(config, ps, 1)'
    elseif config["type"] == "chebyshev"
        start = ChebyshevPolynomial(config, 0)'
        finish = ChebyshevPolynomial(config, 1)'
    end
    A = [start Z Z;
     Z start Z;
     Z Z start;
     finish Z Z;
     Z finish Z;
     Z Z finish]
    b = [xstar; xcurr] 
    return A, b
 end

 

function ClenshawCurtisWeight(config)
    N = config["N"]
    w = zeros(N+1)
    if mod(N,2) == 0
        w[1] = 1/(N^2 - 1)
        w[end] = w[1]
        for s in 1:convert(Int64,N/2)
            wi = 0
            for j in 0:N/2
                if j == 0
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                elseif j == N/2
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                else 
                    wi = wi + 1 /(1-4*j^2) * cos(2*pi*j*s/N);
                end
            end
            w[s+1] = 4/N*wi
            w[N-s+1] = w[s+1]
        end

    else
        w[1] = 1/N^2;
        w[end] = w[1];
        for s = 1:convert(Int64,(N-1)/2)
            wi = 0;
            for j = 0:(N-1)/2
                if j == 0
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                elseif j == N/2
                    wi = wi + 0.5 /(1-4*j^2) * cos(2*pi*j*s/N);
                else 
                    wi = wi + 1 /(1-4*j^2) * cos(2*pi*j*s/N);
                end
            end
            w[s+1] = 4/N*wi;
            w[N-s+1] = 4/N*wi;
        end
    end
    w = w/2
end

evalMetric(config, x) = config["W"][1] + config["W"][2]*x + config["W"][3]*x^2


function computeEnergy(config, ps, C, L, dL) 
    C_reshape = reshape(C, (config["deg"]+1, config["dim"]))'
    X = C_reshape*L
    dX = C_reshape*dL
    # summing at the collcation points
    sum([dX[:,i]'*(evalMetric(config, X[1,i])\dX[:,i])*ps["weights"][i] for i in 1:config["N"]+1])
end


function computeJacobian(config, ps, C, L, dL)
    N = config["N"]
    deg = config["deg"]
    C_reshape = reshape(C, (deg+1, config["dim"]))'
    X = C_reshape*L
    dX = C_reshape*dL
    gx = zeros(deg+1,N+1)
    gy = zeros(deg+1,N+1)
    gz = zeros(deg+1,N+1)
    for j in 1:N+1
        me = X[1,j]
        Wx = evalMetric(config, me)
        M = inv(Wx)
        dW = config["W"][2] + 2*config["W"][3]*me   # differential W wrt x
        Mdot = -M*dW*M
        d = dX[:,j]   # delta
        dsj = ps["weights"][j]   # clencurt weight
        dLds = dL[:,j]*dsj       
        Lds = L[:,j]*dsj
        Md = 2*(M*d)
        gx[:,j] = Md[1]*dLds + (d'*Mdot*d)[1]*Lds;
        gy[:,j] = Md[2]*dLds
        gz[:,j] = Md[3]*dLds;
    end
    g = [gx; gy; gz]
    sparsify(sum(g, 2), config["sparse_eps"])           
end


# dLag = hcat([DiffLagrangePolynomial(config, ps, ni) for ni in ps["nodes"]]...)
# Lag = hcat([LagrangePolynomial(config, ps, ni) for ni in ps["nodes"]]...)

function computeError(config, ps_result)
    orig_N = config["N"]
    config["N"] = 100
    nodes = CGLnodes(config)
    weight = ClenshawCurtisWeight(config)
    if config["type"] == "lagrange"
        dL = hcat([DiffLagrangePolynomial(config, ps, ni) for ni in nodes]...)
        L = hcat([LagrangePolynomial(config, ps, ni) for ni in nodes]...) 
    elseif config["type"] == "chebyshev"
        dL = hcat([DiffChebyshevPolynomial(config, ni) for ni in nodes]...)
        L = hcat([ChebyshevPolynomial(config, ni) for ni in nodes]...) 
    end
    C_reshape = reshape(ps_result["C"], (config["deg"]+1, config["dim"]))'
    X = C_reshape*L
    dX = C_reshape*dL

    V = hcat([dX[:,i]'''*inv(evalMetric(config, X[1,i]))*dX[:,i]'' for i in 1:config["N"]+1]...)
    error = sum((V[i] - ps_result["E0"]).^2*weight[i] for i in 1:config["N"]+1)
    config["N"] = orig_N
    return error
end


# running the optimisation loop
function pseudospectral_geodesic(config, xstar, xcurr)
    ps = ps_params(config)
    C = zeros(config["dim"]*(config["deg"]+1))
    E0 = 0
    tic()
    L = 0
    dL = 0
    for rep in 1:config["repetition"]
        C = zeros(config["dim"]*(config["deg"]+1))
        if config["type"] == "lagrange"
            dL = hcat([DiffLagrangePolynomial(config, ps, ni) for ni in ps["nodes"]]...)
            L = hcat([LagrangePolynomial(config, ps, ni) for ni in ps["nodes"]]...) 
            for q = 1:config["dim"]
                C[[1, config["deg"] + 1, 2*config["deg"] + 3]] = xstar
                C[[config["deg"]+1, 2*(config["deg"] + 1), 3*(config["deg"] + 1)]] = xcurr
            end
        elseif config["type"] == "chebyshev"
            dL = hcat([DiffChebyshevPolynomial(config, ni) for ni in ps["nodes"]]...)
            L = hcat([ChebyshevPolynomial(config, ni) for ni in ps["nodes"]]...) 
            for q = 1:config["dim"]
                C[(config["deg"]+1)*(q-1)+1] = 0.5 * (xstar[q] + xcurr[q]);
                C[(config["deg"]+1)*(q-1)+2] = 0.5 * (xcurr[q] - xstar[q]);
            end
        end
        # initialisation
        A, b = Constraints(config, ps, xstar, xcurr)
        E0 = computeEnergy(config, ps, C, L, dL)
        H = eye(config["dim"]*(config["deg"] + 1))
        H_default = H
        g = computeJacobian(config, ps, C, L, dL)
        # to enter while loop
        step_dir = 1
        alpha = 1   

        while norm(alpha*step_dir) > config["rel_tol"]
            if det(H) < 0
                H = H_default
            end
            KKT_mat = [H A';A zeros(2*config["dim"],2*config["dim"])]
            KKT_sol = KKT_mat\[g ; A*C - b]
            step_dir = -KKT_sol[1:config["dim"]*(config["deg"] + 1)]
            # backtracking line search
            m = g'*step_dir
            t = - config["c"]*m
            alpha = config["alpha0"]
            EE = 0
            while true
                EE = computeEnergy(config, ps, C+alpha*step_dir, L, dL)
                if  (E0 - EE)[1] >= (alpha*t)[1]
                    break
                end
                alpha = config["tau"]*alpha
            end
            E0 = EE
            s = alpha*step_dir  # column
            C = C + s
            g0 = g # column
            g = computeJacobian(config, ps, C, L, dL)
            y = (g - g0) # column
            H = H - (H*(s*s')*H)/(s'*H*s) + (y*y')/(y'*s)
        end
    end
    comp_time = toc()/config["repetition"]
    
    C_reshape = reshape(C, (config["deg"]+1, config["dim"]))'
    X = C_reshape*L
    dX = C_reshape*dL
    V = hcat([dX[:,i]'''*inv(evalMetric(config, X[1,i]))*dX[:,i]'' for i in 1:config["N"]+1]...)
    error = sqrt(sum((V[i] - E0).^2*ps["weights"][i] for i in 1:config["N"]+1))/E0
    
    ps_result = Dict("C" => C,
                     "V" => V,
                     "E0" => E0,
                     "avg_time" => comp_time,
                     "error" => error, 
                     "X" => X,
                     "dX" => dX,
                     "ps" => ps)
    return ps_result


end


    
function computeController(config, xstar, xcurr)
    ps_result = pseudospectral_geodesic(config, xstar, xcurr)
    x1 = ps_result["X"][1,:]
    dX = ps_result["dX"]
    rho = 2 + 1.28210772972104 * x1 + 8.34575059042599 * x1.^2;
    weights = ps_result["ps"]["weights"]
    B = [0,0,1]
    u = -0.5*sum(rho[i]*B'/evalMetric(config, x1[i])*dX[:,i]*weights[i] for i in 1:config["N"]+1)
end

dynamicSystem(x, u) = [-x[1] + x[3] ; x[1]^2 - x[2] - 2 * x[1] * x[3] + x[3] ; -x[2] + u]



function RKIntegrate_CCM(config, xstar, xinit, t_max = 10, dt = 0.05)
    t = 0:dt:t_max
    state = []
    push!(state, xinit)
    xcurr = xinit
    time = []
    u = []
    for i in 1:length(t)-1
        tic()
        u1 = computeController(config, xstar, xcurr)
        dx1 = dynamicSystem(xcurr, u1)
        k1 = dt*dx1

        u2 = computeController(config, xstar, xcurr+0.5*k1)
        dx2 = dynamicSystem(xcurr+0.5*k1, u2)
        k2 = dt*dx2

        u3 = computeController(config, xstar, xcurr+0.5*k2)
        dx3 = dynamicSystem(xcurr+0.5*k2, u3)
        k3 = dt*dx3

        u4 = computeController(config, xstar, xcurr+k3)
        dx4 = dynamicSystem(xcurr+k3, u4)
        k4 = dt*dx4
        xnext =  xcurr + (k1 + 2*k2 + 2*k3 + k4)/6
        push!(state, xnext)
        push!(u, (u1 + 2*u2 + 2*u3 + u4)/6)
        push!(time, toc());
        xcurr = xnext
    end
    state = hcat(state...)'
    u = hcat(u...)'

    return t, state, u, time
end


        

function EulerIntegrate_CCM(config, xstar, xinit, t_max = 10, dt = 0.05)
    t = 0:dt:t_max
    state = []
    push!(state, xinit)
    xcurr = xinit
    comp_time = []
    u = []
    for i in 1:length(t)-1
        tic()
        ui = computeController(config, xstar, xcurr)
        dx = dynamicSystem(xcurr, ui)
        xnext =  xcurr + dx*dt
        push!(state, xnext)
        push!(u, ui)
        push!(comp_time, toc());
        xcurr = xnext
    end
    state = hcat(state...)'
    u = hcat(u...)'
    return t, state, u, comp_time
end







        
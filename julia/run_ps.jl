
include("chebyshev_pseudospectral.jl")


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
              0 0 0]; 

# setting up interpolating nodes CGL


config = Dict("N" => 20,
              "deg" => 7,
              "alpha0" => 1,          # initial step size
              "c" => 0.1,             # termination condition for backtracking line search
              "tau" => 0.1,           # rescale factor for backtracking line search
              "sparse_eps" => 1E-20,
              "rel_tol" => 1E-7,
              "W" => [WC, WL, WQ],
              "repetition" => 100,
              "dim" => 3,
              "type" => "chebyshev"
               )
if config["type"] == "lagrange"
  config["deg"] = config["N"]
end

# finding optimal K_add node
xstar = [0,0,0]
xcurr = [1, 1, 1]*9
cheby_k_d = Dict()
for k in 1:4
  for d in 1:30 
    config["deg"] = d
    config["N"] = config["deg"] + k
    ps_result = pseudospectral_geodesic(config, xstar, xcurr)
    cheby_k_d["k_$(k)_d_$(d)"] = ps_result
    end
end

save("results/cheby_k_add_deg.jld", "cheby_k_d", cheby_k_d)

# for optimal K_add, results for different curr points
curr = collect(1:2:9)
xstar = [0,0,0]

cheby_curr_k = Dict()
for k in curr
  for d in 1:10
    xcurr = [1, 1, 1]*k
    config["deg"] = d
    config["N"] = config["deg"] + 4
    ps_result = pseudospectral_geodesic(config, xstar, xcurr)
    cheby_curr_k["curr_$(k)_k4_d_$(d)"] = ps_result
    end
end

save("results/cheby_curr_fixedk_deg.jld", "cheby_curr_fixedk", cheby_curr_k)

# Chebyshev and Lagrangian polynomial comparison to show we want to use Chebyshev

config = Dict("N" => 20,
              "deg" => 7,
              "alpha0" => 1,          # initial step size
              "c" => 0.1,             # termination condition for backtracking line search
              "tau" => 0.1,           # rescale factor for backtracking line search
              "sparse_eps" => 1E-20,
              "rel_tol" => 1E-7,
              "W" => [WC, WL, WQ],
              "repetition" => 10,
              "dim" => 3,
              "type" => "lagrange"
               )
if config["type"] == "lagrange"
  config["deg"] = config["N"]
end

xstar = [0,0,0]
lagrange_curr = Dict()
for k in 1:6
  for d in 1:10
    xcurr = [1, 1, 1]*k
    config["deg"] = d
    config["N"] = config["deg"]
    ps_result = pseudospectral_geodesic(config, xstar, xcurr)
    lagrange_curr["curr_$(k)_d_$(d)"] = ps_result
    end
end

save("results/lagrange_curr_deg.jld", "lagrange_curr", lagrange_curr)


# cheby deg = N, different curr with different deg
config = Dict("N" => 20,
              "deg" => 7,
              "alpha0" => 1,          # initial step size
              "c" => 0.1,             # termination condition for backtracking line search
              "tau" => 0.1,           # rescale factor for backtracking line search
              "sparse_eps" => 1E-20,
              "rel_tol" => 1E-7,
              "W" => [WC, WL, WQ],
              "repetition" => 10,
              "dim" => 3,
              "type" => "chebyshev"
               )
if config["type"] == "lagrange"
  config["deg"] = config["N"]
end

xstar = [0,0,0]
chebyshev_curr = Dict()
for k in 1:6
  for d in 1:10
    xcurr = [1, 1, 1]*k
    config["deg"] = d
    config["N"] = config["deg"]
    ps_result = pseudospectral_geodesic(config, xstar, xcurr)
    chebyshev_curr["curr_$(k)_d_$(d)"] = ps_result
    end
end

save("results/cheby_curr_deg.jld", "chebyshev_curr", chebyshev_curr)


# complexity numerical simulation


# Plotting and Analysis

# Chebyshev and Lagrangian comparison

cheby_curr_deg = load("results/cheby_curr_deg.jld")
cheby_curr_deg = cheby_curr_deg["chebyshev_curr"]
lagrange_curr_deg = load("results/lagrange_curr_deg.jld")
lagrange_curr_deg = lagrange_curr_deg["lagrange_curr"]

cheby = []
for k in 1:6
    for d in 1:10
        if cheby_curr_deg["curr_$(k)_d_$(d)"]["error"][1] < 1E-5
            er = cheby_curr_deg["curr_$(k)_d_$(d)"]["error"][1]
            time = cheby_curr_deg["curr_$(k)_d_$(d)"]["avg_time"]
            push!(cheby, [er, time, d, k])
            break
        end
    end
end
cheby = hcat(cheby...)'

lagrange = []
for k in 1:6
    for d in 1:10
        if lagrange_curr_deg["curr_$(k)_d_$(d)"]["error"][1] < 1E-5
            er = lagrange_curr_deg["curr_$(k)_d_$(d)"]["error"][1]
            time = lagrange_curr_deg["curr_$(k)_d_$(d)"]["avg_time"]
            push!(lagrange, [er, time, d, k])
            break
        end
    end
end
lagrange = hcat(lagrange...)'

fig = figure(figsize = [7, 4])
plot(cheby[:,4], cheby[:,2], lagrange[:,4], lagrange[:,2], marker = "o", markersize = 10)
ax = gca()
ax[:set_yscale]("log")
legend(["Chebyshev", "Lagrange"], loc = "lower right")
xlabel("Starting Factor", fontsize = 20)
ylabel("Comp. Time (s)", fontsize = 20)
axis([0.5, 3.5, 0.001, 0.1])






# how to choose k_add

cheby_k_d = load("results/cheby_k_add_deg.jld")
cheby_k_d = cheby_k_d["cheby_k_d"]

error = []
time = []
for d in 1:30 
    error_d = []
    time_d = []
    for k in 2:4
        push!(error_d, cheby_k_d["k_$(k)_d_$(d)"]["error"][1])
        push!(time_d, cheby_k_d["k_$(k)_d_$(d)"]["avg_time"][1])
    end
    push!(error, error_d)
    push!(time, time_d)
end
error = hcat(error...)
time = hcat(time...)

fig = figure(figsize = (16,3))
subplot(121)
plot(error', marker = "o", markersize = 8)
ax = gca()
ax[:set_yscale]("log")
# legend(["a=2", "a=3", "a=4"], fontsize = 20)
xlabel("Degree (D)", fontsize = 20)
ylabel("\$\\mathcal{E}\$", fontsize = 20)
subplot(122)
plot(time', marker = "o")
ax = gca()
ax[:set_yscale]("log")
legend(["a=2", "a=3", "a=4"], fontsize = 20, loc = "lower right")
xlabel("Degree (D)", fontsize = 20)
ylabel("Comp. Time (s)", fontsize = 20)



curr = collect(1:2:9)
xstar = [0,0,0]

cheby_curr_k = Dict()
for k in curr
  for d in 1:10
    xcurr = [1, 1, 1]*k
    config["deg"] = d
    config["N"] = config["deg"] + 4
    ps_result = pseudospectral_geodesic(config, xstar, xcurr)
    cheby_curr_k["curr_$(k)_k4_d_$(d)"] = ps_result
    end
end

save("results/cheby_curr_fixedk_deg.jld", "cheby_curr_fixedk", cheby_curr_k)



#  for optimal K_add
curr = collect(1:2:9)
cheby_curr = load("results/cheby_curr_fixedk_deg.jld")
data = []
for k in curr
    for d in 1:10
        if cheby_curr["cheby_curr_fixedk"]["curr_$(k)_k4_d_$(d)"]["error"][1] < 1E-6
            er = cheby_curr["cheby_curr_fixedk"]["curr_$(k)_k4_d_$(d)"]["error"][1]
            ti = cheby_curr["cheby_curr_fixedk"]["curr_$(k)_k4_d_$(d)"]["avg_time"]
            push!(data, [k, d, er, ti])
            break
        end
    end
end
data = hcat(data...)'


# 5×4 Array{Float64,2}:
#  1.0  4.0  1.82729e-8  0.00811889
#  3.0  4.0  2.4745e-7   0.00848628
#  5.0  5.0  3.96295e-7  0.0148703 
#  7.0  6.0  1.91516e-7  0.0206947 
#  9.0  7.0  9.09611e-8  0.0254276 


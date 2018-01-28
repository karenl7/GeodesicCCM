function LQRSimulate(xinit, t_max = 10, dt = 0.05)
    A = [-1 0 1 ; 0 -1 1 ; 0 -1 0];
    B = [0;0;1];
    Q = eye(3);
    R = 1;
    K = [0.2000   -0.2000    1.0000]

    dx(x, u) =  [-x[1] + x[3] ; x[1]^2-x[2]-2*x[1]*x[3]+x[3] ; -x[2] + u] 

    t = 0:dt:t_max
    state = []
    push!(state, xinit)
    xcurr = xinit
    time = []
    u = []
    for i in 1:length(t)-1
    tic()
    u1 = -K*xcurr
    dx1 = dx(xcurr, u1)
    k1 = dt*dx1

    u2 = -K*(xcurr+0.5*k1)
    dx2 = dx(xcurr+0.5*k1, u2)
    k2 = dt*dx2

    u3 = -K*(xcurr+0.5*k2)
    dx3 = dx(xcurr+0.5*k2, u3)
    k3 = dt*dx3

    u4 = -K*(xcurr+k3)
    dx4 = dx(xcurr+k3, u4)
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

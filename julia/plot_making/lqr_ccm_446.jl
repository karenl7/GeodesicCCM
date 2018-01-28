xinit = [4, 4, 6]
xstar = [0, 0, 0]
t_lqr, state_lqr, u_lqr, time_lqr = LQRSimulate(xinit, 10);
t_ccm, state_ccm, u_ccm, time_ccm = EulerIntegrate_CCM(config, xstar, xinit, 10, 0.01);
t_ccm_9, state_ccm_9, u_ccm_9, time_ccm_9 = EulerIntegrate_CCM(config, xstar, [9, 9, 9], 10, 0.01);




fig = figure(figsize = (7,3.5))
plot(t_lqr'', state_lqr[:,1], color = "g")
plot(t_ccm'', state_ccm[:,1], color = "r")
plot(t_ccm_9'', state_ccm_9[:,1], color = "b")
plot(t_lqr'', state_lqr[:,2:3], color = "g")
plot(t_ccm'', state_ccm[:,2:3], color = "r")
plot(t_ccm_9'', state_ccm_9[:,2:3], color = "b")

axis([0, 10, -60, 20])
legend(["LQR", "CCM", "CCM_9"], loc = "lower right")
xlabel("Sim. Time (s)", fontsize = 20)
ylabel("\$x_1, x_2, x_3\$", fontsize = 20)


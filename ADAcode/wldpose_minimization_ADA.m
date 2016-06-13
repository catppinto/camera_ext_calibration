function F = wldpose_minimization_ADA(x, wld_pose, cam_pose, x0)
alpha = 1;
load ADA_Data;

k= T_e2c(1:3, 1:3);
t_Est = [x(1); x(2) ; x(3)];
k = [k t_Est; 0 0 0 1];

wld_estimate = k*cam_pose; 

error = sqrt((wld_pose(1) - wld_estimate(1))^2 + ...
    (wld_pose(2) - wld_estimate(2))^2 + ...
    (wld_pose(3) - wld_estimate(3))^2 );
d = sqrt((t_Est(1) - x0(1))^2 + ...
    (t_Est(2) - x0(2))^2 + ...
    (t_Est(3) - x0(3))^2 );
F = 10*error + alpha * d;

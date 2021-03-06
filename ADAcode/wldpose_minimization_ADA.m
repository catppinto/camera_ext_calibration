function F = wldpose_minimization_ADA(x, wld_pose, cam_pose, K_ec, transforms)

T_rb2e = transforms.T_rb2e;
T_w2rb = transforms.T_w2rb;

%%
alpha = [ 0.01 0.01 0.1]';
beta = [0.1 0.1 0.001]';

k= K_ec(1:3, 1:3);
t_Est = [x(1); x(2); x(3)];
k = [k t_Est; 0 0 0 1];

t_w = T_w2rb * T_rb2e * k*cam_pose; 

if(size(wld_pose, 2) ==1 ) 
    error = sqrt((wld_pose(1) - t_w(1))^2 + ...
        (wld_pose(2) - t_w(2))^2 + ...
        (wld_pose(3) - t_w(3))^2 );
    d = sqrt((t_Est(1) - K_ec(1,4))^2 + ...
        (t_Est(2) - K_ec(2,4))^2  + ...
        (t_Est(3) - K_ec(3,4))^2 );
    F = error + alpha * d;
else
    error = [abs( (wld_pose(1,:) - t_w(1, :))) ; 
             abs( (wld_pose(2, :) - t_w(2,:))) ;
             abs( (wld_pose(3, :) - t_w(3, :)))];
    error = sum(error,2);
    d = [abs(t_Est(1) - K_ec(1,4)).^2; 
         abs(t_Est(2) - K_ec(2,4)).^2;
         abs(t_Est(3) - K_ec(3,4)).^2];
    F = sum(beta.*error) + sum(alpha .* d);
    
end

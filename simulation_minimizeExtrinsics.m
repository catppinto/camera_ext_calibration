function [F] = simulation_minimizeExtrinsics(x, tagposes,  corners, K_ext0)

zhat_table = [0 0 1]';
z_table = 0.5;   

alpha = 10;
delta = 10;

A = [529.2945         0  466.9604         0; ...
            0  531.2835  273.2594         0; ...
            0         0    1.0000         0; ...
            0         0         0    1.0000];

%% Conversion to R,t extrinsiccs 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

K_wc_est = [R_wc t_wc ; 0 0 0 1];


%% Load Tag Poses

counter = 0;
for i=1:size(tagposes,3) 
        
    t_c = tagposes(:, :, i);
    t_w = K_wc_est * t_c;
    t_w_estimates(:, :, i) = t_w;
    
    zhat_tag = t_w(1:3, 3);
    
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_z(counter) = norm(z_tag - z_table);
    
    CosTheta = dot(zhat_tag,zhat_table)/(norm(zhat_tag)*norm(zhat_table));
    theta_error(counter) = acos(CosTheta);
       
end

% norm extrinsics
gamma = [0.0001 10 10 0.01]';
t_ext0 = K_ext0(4:6); 
norm_tExt = [abs(t_ext0(1) - x(4)).^2; 
     abs(t_ext0(2) - x(5)).^2;
     abs(t_ext0(3) - x(6)).^2];
r_ext0 = K_ext0(1:3); 
norm_rExt = norm(r_ext0 - x(1:3));
norm_ext = [norm_rExt ; norm_tExt];

%% corners

tag_size = 0.06;
ts = tag_size/2;
for i=1:size(t_w_estimates,3)
    t = t_w_estimates(:, :, i);
    corner1_est = t; corner1_est(1,4) = t(1,4)+ts; corner1_est(2,4) = t(2,4)+ts; 
    corner2_est = t; corner2_est(1,4) = t(1,4)-ts; corner2_est(2,4) = t(2,4)+ts;
    corner3_est = t; corner3_est(1,4) = t(1,4)-ts; corner3_est(2,4) = t(2,4)-ts;
    corner4_est = t; corner4_est(1,4) = t(1,4)+ts; corner4_est(2,4) = t(2,4)-ts;
    
    T_w1(:, :, i) = corner1_est; 
    T_w2(:, :, i) = corner2_est;
    T_w3(:, :, i) = corner3_est;
    T_w4(:, :, i) = corner4_est;
end

for n =1:size(t_w_estimates,3)
   T_c1(:, :, n) = inv(K_wc_est) * T_w1(:, :, n);
   T_c2(:, :, n) = inv(K_wc_est) * T_w2(:, :, n);
   T_c3(:, :, n) = inv(K_wc_est) * T_w3(:, :, n);
   T_c4(:, :, n) = inv(K_wc_est) * T_w4(:, :, n);
end
        
c1_est = []; c2_est= []; c3_est= []; c4_est= [];
for n =1:size(t_w_estimates,3)
   tc_1 = T_c1(:, :, n) ;
   corner1_est = A * tc_1(1:4, 4);
   c1_est(n, :) = [corner1_est(1)/corner1_est(3) corner1_est(2)/corner1_est(3) 1];

   tc_2 = T_c2(:, :, n) ;
   corner2_est = A * tc_2(1:4, 4);
   c2_est(n, :) = [corner2_est(1)/corner2_est(3) corner2_est(2)/corner2_est(3) 1];
   
   tc_3 = T_c3(:, :, n) ;
   corner3_est = A * tc_3(1:4, 4);
   c3_est(n, :) = [corner3_est(1)/corner3_est(3) corner3_est(2)/corner3_est(3) 1];
   
   tc_4 = T_c4(:, :, n) ;
   corner4_est = A * tc_4(1:4, 4);
   c4_est(n, :) = [corner4_est(1)/corner4_est(3) corner4_est(2)/corner4_est(3) 1];
end

for n =1:size(t_w_estimates,3)
    norm_c1(n) = norm(c1_est(n, :) - corners.c1(n, :));
    norm_c2(n) = norm(c2_est(n, :) - corners.c2(n, :));
    norm_c3(n) = norm(c3_est(n, :) - corners.c3(n, :));
    norm_c4(n) = norm(c4_est(n, :) - corners.c4(n, :));
end

reprojection_error = sum(norm_c1) + sum(norm_c2) + sum(norm_c3) + sum(norm_c4);


%% cost function 
F = 5*sum(norm_z) + sum(alpha*theta_error) + 10*sum(gamma.*norm_ext) + delta*reprojection_error; 




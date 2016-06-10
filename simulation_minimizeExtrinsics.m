function [F] = simulation_minimizeExtrinsics(x, tagposes)

zhat_table = [0 0 1]';
z_table = 0.5;   

alpha = 0.03;

%% Conversion to R,t extrinsiccs 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

K_wc_est = [R_wc t_wc ; 0 0 0 1];


%% Load Tag Poses

counter = 0;
for i=1:size(tagposes,3) 
        
    t_c = tagposes(:, :, i);
    t_w = K_wc_est * t_c;
    
    axis_angle = vrrotmat2vec(t_w(1:3, 1:3));
    if(axis_angle(3) <0)
        zhat_tag = -1 .* axis_angle(1:3);
    else
        zhat_tag = axis_angle(1:3);
    end
    
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_z(counter) = norm(z_tag - z_table);
    norm_zhat(counter) = dot(zhat_tag,zhat_table);
       
end
    
F = sum(norm_z) - sum(alpha*norm_zhat);
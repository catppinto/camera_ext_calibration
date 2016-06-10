function [F] = minimize_tableTagErrorADA(x)

%% get projection data from HERB (init estimate) 
load ADA_Data; %% A, K, P

zhat_table = [0 0 1]';
z_table = tag_height;   

alpha = 0.003;

%% Conversion to R,t extrinsiccs 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

K_wc_est = [R_wc t_wc ; 0 0 0 1];


%% Load Tag Poses

load('/home/cat/Documents/CMU_Herb/camera_ext_calibration/dataFromADA/ADAtags_08062016.mat')
% only use tags_test

for i=1:size(tags_train,1) 
    
    t = tags_train(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    % for ADA
    t_w = T_tb2rb*T_rb2e * K_wc_est * t_c;
    
    zhat_tag = t_w(1:3, 3);
    z_tag =  t_w(3, 4);
        
    norm_z(i) = norm(z_tag - z_table);
    norm_zhat(i) = dot(zhat_tag,zhat_table);
       
end

    
F = sum(norm_z) + sum(alpha*norm_zhat);


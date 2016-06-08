function [F] = minimize_tableTagError(x)

%% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P
load HERB_tableData.mat

zhat_table = [0 0 1]';
z_table = tb_height;   

alpha = 0.2;

%% Conversion to R,t extrinsiccs 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

K_wc_est = [R_wc t_wc ; 0 0 0 1];


%% Load Tag Poses

load('dataFromTableRotationTest/tagTableRotationTest07062016.mat')
% only use tags_test

for i=1:size(tags_train,1) 
    
    
    t = tags_train(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    marker_rot_tabl = (tb_kinbody_offset(1:3, 1:3));
    t_c_ = [inv(marker_rot_tabl)*t_c(1:3, 1:3) t_c(1:3,4) ; 0 0 0 1];
    t_w = K_wc_est * t_c_;
    
    zhat_tag = t_w(1:3, 3);
    z_tag =  t_w(3, 4);
    
        
    norm_z(i) = norm(z_tag - z_table);
    norm_zhat(i) = norm(dot(zhat_tag,zhat_table));
    
    
end

    
F = sum(norm_z) - sum(alpha*norm_zhat);


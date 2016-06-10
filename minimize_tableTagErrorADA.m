function [F] = minimize_tableTagErrorADA(x)

%% get projection data from HERB (init estimate) 
load ADA_Data; %% A, K, P

zhat_table = [0 0 1]';
z_table = tag_height;   

alpha = 10;

%% Conversion to R,t extrinsiccs 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

K_wc_est = [R_wc t_wc ; 0 0 0 1];


%% Load Tag Poses

load('dataFromADA/ADAtags_08062016.mat')
% only use tags_test
counter = 0;
for i=1:size(tags_train,1) 
    
    t = tags_train(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    % for ADA
    t_w = T_tb2rb*T_rb2e * K_wc_est * t_c;
    qw = rotationMatrixToQ(t_w(1:3, 1:3))';
    
    axis_angle = vrrotmat2vec(t_w(1:3, 1:3));
    if(axis_angle(3) <0)
        zhat_tag = -1 .* axis_angle(1:3);
    else
        zhat_tag = axis_angle(1:3);
    end
    
    if(zhat_tag(3) > 0.5)
        z_tag =  t_w(3, 4);
        counter = counter + 1;
        norm_z(counter) = norm(z_tag - z_table);
        norm_zhat(counter) = dot(zhat_tag,zhat_table);
    end
       
end
    
F = sum(norm_z) - sum(alpha*norm_zhat);


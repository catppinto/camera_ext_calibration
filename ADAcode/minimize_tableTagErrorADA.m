function [F] = minimize_tableTagErrorADA(x, tagsPose, K_ext0)

%% get projection data from HERB (init estimate) 
load ADA_Data; %% A, K, P
load T_rb2e_11062016;
T_rb2e = T_rb2e_11062016;

zhat_table = [0 0 1]';
z_table = tag_height;   

alpha = 10;
gamma = 0.05;
beta = 0.05;

%% Conversion to R,t extrinsiccs 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

K_wc_est = [R_wc t_wc ; 0 0 0 1];


%% Load Tag Poses

% only use tags_test
counter = 0;
for i=1:size(tagsPose,1) 
    
    t = tagsPose(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    % for ADA
    t_w = T_tb2rb*T_rb2e * K_wc_est * t_c;
    
    zhat_tag = t_w(1:3, 3);
    
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_z(counter) = norm(z_tag - z_table);
    
    CosTheta = dot(zhat_tag,zhat_table)/(norm(zhat_tag)*norm(zhat_table));
    theta_error(counter) = acos(CosTheta);
       
end
    
r_ext0 = K_ext0(1:3); 
norm_rExt = norm(r_ext0 - x(1:3));
t_ext0 = K_ext0(4:6); 
norm_tExt = norm(t_ext0 - x(4:6));
%% cost function 
F = sum(norm_z) + sum(alpha*theta_error)+ gamma*norm_rExt + beta*norm_tExt; 


 
%%%% TODO 
% add Intrinsics from ADA
% add corners (redo calibration)
% Confirm computation of Tw
%%%%

clear
clc
 
%% TRAINING

%% LOADING DATA
load ADA_Data; %% A, K, P
load('dataFromADA/ADAtags_11062016.mat')
corners = [];

%%  find x

K_wc = T_e2c;

w=rodrigues(K_wc(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4)]

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@minimize_tableTagErrorADA,x0,[],[],[],[],[],[],[], options, tags_train, corners, x0);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

K_wc

%% Parameter Refinement

zhat_table = [0 0 1]';
z_table = tag_height; 

% FROM OPTIMIZATION
counter=0;
zhat_tag = [];
for i=1:size(tags_train,1) 
     
    t = tags_train(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    t_w = T_tb2rb*T_rb2e * KK * inv(t_c);
    
    zhat_tag = t_w(1:3, 3);
    
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_zO(counter) = norm(z_tag - z_table);
    
    CosTheta = dot(zhat_tag,zhat_table)/(norm(zhat_tag)*norm(zhat_table));
    theta_errorO(counter) = acos(CosTheta);
   
end

disp('****************************')
disp('Optimized Results')
norm_zO
theta_errorO
disp('****************************')

zhat_table = [0 0 1]';
counter = 0;
for i=1:size(tags_train,1) 
    
    t = tags_train(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    % for ADA
    t_w = T_tb2rb*T_rb2e * K_wc * t_c;
    qw = rotationMatrixToQ(t_w(1:3, 1:3))';
    
    zhat_tag = t_w(1:3, 3);
    
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_z(counter) = norm(z_tag - z_table);
    
    CosTheta = dot(zhat_tag,zhat_table)/(norm(zhat_tag)*norm(zhat_table));
    theta_error(counter) = acos(CosTheta);
 
end

disp('****************************')
disp('Real Results')
norm_zs
theta_error
disp('****************************')


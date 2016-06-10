clear
clc
 
%% TRAINING

addpath(strcat(pwd, '/bkg_code'))
addpath(strcat(pwd, '/dataFromTableRotationTest'))

%% LOADING DATA
load ADA_Data; %% A, K, P
load('/home/cat/Documents/CMU_Herb/camera_ext_calibration/dataFromADA/ADAtags_08062016.mat')


%%  find x

K_wc = T_e2c;

w=rodrigues(K_wc(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4)]

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point', 'StepTolerance', 1e-20);

[x, fval, exittag] = fmincon(@minimize_tableTagErrorADA,x0,[],[],[],[],[],[],[], options);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

K_wc

%% Parameter Refinement

zhat_table = [0 0 1]';
z_table = tag_height; 

% FROM OPTIMIZATION
for i=1:size(tags_train,1) 
     
    t = tags_train(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    t_w = T_tb2rb*T_rb2e * KK * t_c;
        
    zhat_tag = t_w(1:3, 3);
    z_tag =  t_w(3, 4);
        
    norm_z(i) = norm(z_tag - z_table);
    norm_zhat(i) = dot(zhat_tag,zhat_table);
    
end

disp('****************************')
disp('Optimized Results')
norm_z
norm_zhat
disp('****************************')
% %% TESTING
% 
% zhat_table = [0 0 1]';
% z_table = tag_height; 
% 
% % FROM OPTIMIZATION
% for i=1:size(tags_test,1) 
%      
%     t = tags_test(i, :);
%     
%     t_c = eye(4,4); 
%     t_c(1:3, 4) = t(1:3);
%     q = t(4:7);
%     t_c(1:3, 1:3) = quatToRotationMatrix(q);
%     
%     t_w = T_tb2rb*T_rb2e * KK * t_c;
%         
%     zhat_tag = t_w(1:3, 3);
%     z_tag =  t_w(3, 4);
%         
%     norm_z(i) = norm(z_tag - z_table);
%     norm_zhat(i) = dot(zhat_tag,zhat_table);
%     
% end
% 
% disp('****************************')
% disp('Optimized Results')
% norm_z
% norm_zhat
% disp('****************************')
% 
% % FROM INITIAL GUESS
% for i=1:size(tags_test,1) 
%      
%     t = tags_test(i, :);
%     
%     t_c = eye(4,4); 
%     t_c(1:3, 4) = t(1:3);
%     q = t(4:7);
%     t_c(1:3, 1:3) = quatToRotationMatrix(q);
%     
%     t_w = T_tb2rb*T_rb2e * K_wc * t_c;
%         
%     zhat_tag = t_w(1:3, 3);
%     z_tag =  t_w(3, 4);
%         
%     norm_z_0(i) = norm(z_tag - z_table);
%     norm_zhat_0(i) = dot(zhat_tag,zhat_table);
%     
% end
% 
% disp('****************************')
% disp('Initial Guess')
% norm_z_0
% norm_zhat_0
% disp('****************************')
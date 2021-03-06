clear
clc
path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/';
path = '~/cat_workspace/src/ExtrinsicsCalibration/matlab_code/';
addpath([path, 'bkg_code'])
%addpath('/homes/apirespi/cat_workspace/src/ExtrinsicsCalibration/matlab_code/bkg_code')


tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];
%% TRAINING

%% LOADING DATA
load ADA_Data; %% A, K, P
%load([path, 'dataFromADA/ADAtags_13062016_wldPoseTag.mat'])

load([path, 'dataFromADA/ADAtags_14062016_wldPoseTag.mat'])
%%  find x

K_ec = [   0.9988    0.0305   -0.0385   -0.0470;
    0.0390   -0.9686    0.2457    0.0500;
   -0.0298   -0.2469   -0.9686   -0.0891;
         0         0         0    1.0000]
     
w=rodrigues(K_ec(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[K_ec(1,4), K_ec(2,4), K_ec(3,4)];

wld_pose = [tags.worldpose(:, :) ones(size(tags.worldpose,1), 1)]';
cam_pose = [tags.pose(:, 1:3) ones(size(tags.worldpose,1), 1)]';

% % simulation  
% K_ec_goal = [   0.9988    0.0305   -0.0385   0.134;
%                 0.0390   -0.9686    0.2457    0.21;
%                -0.0298   -0.2469   -0.9686   -0.0891;
%                      0         0         0    1.0000]
% 
% cam_pose = inv(T_tb2rb * T_rb2e * K_ec_goal)*wld_pose; 
     


%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@wldpose_minimization_ADA,x0,[],[],[],[],[],[],[], options, wld_pose, cam_pose, K_ec);

% get refined extrinsics
t_est = x(1:3)';

k_ec_opt = [K_ec(1:3, 1:3), t_est ; 0 0 0 1]

K_ec

t_w_opt = T_tb2rb * T_rb2e * k_ec_opt*cam_pose 

t_w_prior_opt = T_tb2rb * T_rb2e * K_ec*cam_pose

wld_pose


%% error of projection 


T_e2c = [    1.0000         0         0   -0.0470;
         0   -0.9755    0.2198    0.0500;
         0   -0.2198   -0.9755   -0.0900;
         0         0         0    1.0000];
t_w_est1 = T_tb2rb * T_rb2e * T_e2c*cam_pose 

[mean_abs_error_i, std_deviation_i,  mean_error_3D_i, std_deviation_3D_i] = translationErrorBetweenPointsInWorld(wld_pose, t_w_est1)


[mean_abs_error_prioropt, std_deviation_prioropt,  mean_error_3D_prioropt, std_deviation_3D_prioropt] = translationErrorBetweenPointsInWorld(wld_pose, t_w_prior_opt)

[mean_abs_error_opt, std_deviation_opt,  mean_error_3D_opt, std_deviation_3D_opt] = translationErrorBetweenPointsInWorld(wld_pose, t_w_opt)






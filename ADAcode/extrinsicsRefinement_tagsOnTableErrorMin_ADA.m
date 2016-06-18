%%%% TODO 
% add Intrinsics from ADA
% add corners (redo calibration)
% Confirm computation of Tw
%%%%

clear
clc
% path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/';
% addpath([path, 'bkg_code'])
path = '~/cat_workspace/src/ExtrinsicsCalibration/matlab_code/';
addpath([path 'bkg_code'])
tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];
%% LOADING DATA
load ADA_Data; %% A, K, P
%load([path, 'dataFromADA/ADAtags_13062016_trueSizeTag.mat'])
load([path, 'dataFromADA/ADA_depthtests.mat'])
load([path, 'dataFromADA/ADA_depthtests_gt.mat'], 'T_rb2e')
tags_train = tags; %% change this 

corners = [];

disp('Perform Physics Based Refinement')
%%  find x

init_K =    [ 1.0000    0.0000    0.0000   -0.3995
    0.0000   -0.9639    0.2663    0.4263
   -0.0000   -0.2663   -0.9639    0.7148
         0         0         0    1.0000];

w=rodrigues(init_K(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, init_K(1,4), init_K(2,4), init_K(3,4)]

%% optimization using fmincon
options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@minimize_tableTagErrorADA,x0,[],[],[],[],[],[],[], options, tags_train.pose, x0);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

PB_K = [R_est, t_est ; 0 0 0 1]

init_K

%% Parameter Refinement

zhat_table = [0 0 1]';
z_table = tag_height; 

A = [567.7969970703125, 0.0, 319.5, 0.0; ...
     0.0, 567.7969970703125, 239.5, 0.0;...
     0.0, 0.0, 1.0, 0.0;...
     0.0, 0.0, 0.0, 1.0];
 
% FROM REAL data
counter=0;
zhat_tag = [];


fullcampose = false;
if(length(size(tags_train.pose)) > 2)
    fullcampose = true;
end
if(fullcampose)
    s = size(tags_train,3); 
else
    s = size(tags_train,1);
end

 if(fullcampose)
    tag_rotation_offset = eye(4);
else
    tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];
 end

for i=1:s
     
           
    if(fullcampose)
        t = tags_train.pose(:,  :, 1);
        t_c = t;
    else
        t = tags_train.pose(i, :);
        t_c = eye(4,4); 
        t_c(1:3, 4) = t(1:3);
        q = t(4:7);
        t_c(1:3, 1:3) = quatToRotationMatrix(q);
    end
    
    t_wr = T_tb2rb * T_rb2e * init_K * t_c;
    t_wr = tag_rotation_offset * t_wr;
    fprintf('real: %3.7f \n', t_wr(3,4)); 
    t_w_estimatesr(:, :, i) =t_wr;
    
    t_wo = T_tb2rb * T_rb2e * PB_K * t_c;
    t_wo = tag_rotation_offset * t_wo;
    fprintf('opt: %3.7f \n\n', t_wo(3,4)); 
    t_w_estimateso(:, :, i) =t_wo;
    
    zhat_tagr = t_wr(1:3, 3);
    z_tagr =  t_wr(3, 4);
    
    zhat_tago = t_wo(1:3, 3);
    z_tago =  t_wo(3, 4);
    
    counter = counter + 1;
    norm_zr(counter) = norm(z_tagr - z_table);
    CosTheta = dot(zhat_tagr,zhat_table)/(norm(zhat_tagr)*norm(zhat_table));
    theta_errorr(counter) = acos(CosTheta);
    
    counter = counter + 1;
    norm_zo(counter) = norm(z_tago - z_table);
    CosTheta = dot(zhat_tago,zhat_table)/(norm(zhat_tago)*norm(zhat_table));
    theta_erroro(counter) = acos(CosTheta);
   
end


%% Perform World Known Poses Refinement
disp('Perform World Known Poses Refinement')
% load([path, 'dataFromADA/ADAtags_14062016_wldPoseTag.mat'])
load([path, 'dataFromADA/ADA_depthtests_gt.mat'])
%%  find x

w=rodrigues(PB_K(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[PB_K(1,4), PB_K(2,4), PB_K(3,4)];

wld_pose = [tags.worldpose(:, :) ones(size(tags.worldpose,1), 1)]';

if(fullcampose)
    for t = 1:size(tags.pose,3)
        tag = tags.pose(:, :, t);
        cam_pose(:, t) = tag(:, 4);
    end
else
    cam_pose = [tags.pose(:, 1:3) ones(size(tags.worldpose,1), 1)]';
end
%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@wldpose_minimization_ADA,x0,[],[],[],[],[],[],[], options, wld_pose, cam_pose, PB_K);

% get refined extrinsics
t_est = x(1:3)';

WK_K = [PB_K(1:3, 1:3), t_est ; 0 0 0 1];

WK_K

%% show results

% matrices
init_K

PB_K

WK_K


% world poses

wld_pose

t_w_init = T_tb2rb * T_rb2e * init_K *cam_pose

t_w_PB = T_tb2rb * T_rb2e * PB_K * cam_pose

t_w_WK = T_tb2rb * T_rb2e * WK_K * cam_pose  


%% errors 

[mean_abs_error_i, std_deviation_i,  mean_error_3D_i, std_deviation_3D_i] = translationErrorBetweenPointsInWorld(wld_pose, t_w_init)

[mean_abs_error_prioropt, std_deviation_prioropt,  mean_error_3D_prioropt, std_deviation_3D_prioropt] = translationErrorBetweenPointsInWorld(wld_pose, t_w_PB)

[mean_abs_error_opt, std_deviation_opt,  mean_error_3D_opt, std_deviation_3D_opt] = translationErrorBetweenPointsInWorld(wld_pose, t_w_WK)









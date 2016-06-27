clear
clc
close all

%% needed variables (code conflicts)

global is_cmu;
is_cmu = true;

global load_depth;
load_depth = true;

%% Set path 

if(is_cmu)
   path = '~/cat_workspace/src/ExtrinsicsCalibration/matlab_code/';
else
   path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/'; 
end
addpath([path 'bkg_code'])

%% LOADING DATA
load ADA_Data; %% A, K, P

if load_depth
    load([path, 'dataFromADA/ADA_depthtests.mat'])
    load([path, 'dataFromADA/ADA_depthtests_gt.mat'], 'T_rb2e')
    tags_train = tags; 
else
    load([path, 'dataFromADA/ADAtags_13062016_trueSizeTag.mat'])
end

%% matrices needed

T_w2rb = [1., 0., 0., 0.409; 0., 1., 0., 0.338;  0., 0., 1., 0.795;  0., 0., 0., 1.];

if(load_depth) 
    tag_height = T_w2rb(3,4)-(0.0370-0.03);
else
    tag_height = T_w2rb(3,4)-(0.0370-0.003);
end

T_rb2tb = inv(T_tb2rb);


A = [567.7969970703125, 0.0, 319.5, 0.0; ...
     0.0, 567.7969970703125, 239.5, 0.0;...
     0.0, 0.0, 1.0, 0.0;...
     0.0, 0.0, 0.0, 1.0];

transforms = struct('T_w2rb', T_w2rb, 'T_rb2e', T_rb2e , 'T_tb2rb', T_tb2rb, ...
                    'A', A, ...
                    'tag_height', tag_height);

%% Start procedure
corners = [];

disp('Perform Physics Based Refinement')

init_K = T_e2c;
        
w=rodrigues(init_K(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, init_K(1,4), init_K(2,4), init_K(3,4)]

%% optimization using fmincon
options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@minimize_tableTagErrorADA,x0,[],[],[],[],[],[],[], options, tags_train.pose, x0, transforms);

%% 
% get refined extrinsics
R_est = rodrigues(x(1:3));
R_est(3,2) = - R_est(3,2);
R_est(2,3) = - R_est(2,3);

t_est =[x(4);x(5);x(6)];

PB_K = [R_est, t_est ; 0 0 0 1]

init_K

%% Parameter Refinement

zhat_table = [0 0 -1]';
z_table = tag_height; 
 
% FROM REAL data
counter=0;
zhat_tag = [];


fullcampose = false;
if(load_depth)
    fullcampose = true;
end
if(load_depth)
    s = size(tags_train.pose,3); 
else
    s = size(tags_train.pose,1);
end

for i=1:s
             
    if(fullcampose)
        t = tags_train.pose(:,  :, i);
        t_c = t;
    else
        t = tags_train.pose(i, :);
        t_c = eye(4,4); 
        t_c(1:3, 4) = t(1:3);
        q = t(4:7);
        t_c(1:3, 1:3) = quatToRotationMatrix(q);
    end
    
    t_wr = T_w2rb * T_rb2e * init_K * t_c;
    t_w_estimatesr(:, :, i) =t_wr;
    
    t_wo = T_w2rb * T_rb2e * PB_K * t_c;
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

zhat_tagr
zhat_tago

%% Perform World Known Poses Refinement
disp('Perform World Known Poses Refinement')

if(load_depth)
    load([path, 'dataFromADA/ADA_depthtests_gt.mat'])
else
    load([path, 'dataFromADA/ADAtags_14062016_wldPoseTag.mat'])
end

%%  find x

wld_pose = [tags.worldpose(:, :) ones(size(tags.worldpose,1), 1)]';

if(fullcampose)
    for t = 1:size(tags.pose,3)
        tag = tags.pose(:, :, t);
        cam_pose(:, t) = tag(:, 4);
    end
else
    cam_pose = [tags.pose(:, 1:3) ones(size(tags.worldpose,1), 1)]';
end


tag_in_world_gt = T_w2rb * T_rb2tb * wld_pose ;


%% show results

% matrices
init_K

PB_K


% world poses

tag_in_world_gt

t_w_init =  T_w2rb * T_rb2e * init_K * cam_pose

t_w_PB =  T_w2rb * T_rb2e * PB_K * cam_pose


%% errors 

[mean_abs_error_i, std_deviation_i,  mean_error_3D_i, std_deviation_3D_i] = ...
    translationErrorBetweenPointsInWorld(tag_in_world_gt, t_w_init);

[mean_abs_error_PB, std_deviation_PB,  mean_error_3D_PB, std_deviation_3D_PB] = ...
    translationErrorBetweenPointsInWorld(tag_in_world_gt, t_w_PB);

rodrigues(PB_K(1:3 , 1:3))

rodrigues(init_K(1:3 , 1:3))





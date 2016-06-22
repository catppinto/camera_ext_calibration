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

tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];

A = [567.7969970703125, 0.0, 319.5, 0.0; ...
     0.0, 567.7969970703125, 239.5, 0.0;...
     0.0, 0.0, 1.0, 0.0;...
     0.0, 0.0, 0.0, 1.0];

T_w2e =  T_rb2e;
clear T_rb2e ;

transforms = struct('T_w2rb', T_w2rb, 'T_w2e', T_w2e , 'T_tb2rb', T_tb2rb, ...
                    'tag_rotation_offset', tag_rotation_offset, 'A', A, ...
                    'tag_height', tag_height);

%% Start procedure
corners = [];

disp('Perform Physics Based Refinement')

% recover x
init_K = T_e2c;

init_K =    [ 1.0000    0.0000    0.0000   -0.3995;...
    0.0000   -0.9639    0.2663    0.4263;...
   -0.0000   -0.2663   -0.9639    0.7148;...
         0         0         0    1.0000];
         
w=rodrigues(init_K(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, init_K(1,4), init_K(2,4), init_K(3,4)]

%% optimization using fmincon
options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@minimize_tableTagErrorADA,x0,[],[],[],[],[],[],[], options, tags_train.pose, x0, transforms);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

PB_K = [R_est, t_est ; 0 0 0 1]

init_K

%% Parameter Refinement

zhat_table = [0 0 1]';
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

if(fullcampose)
    tag_rotation_offset = eye(4);
else
    tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];
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
    
    t_wr = T_w2e * init_K * t_c;
    t_wr = [tag_rotation_offset(1:3,1:3) * t_wr(1:3,1:3) t_wr(1:3,4); 0 0 0 1];
    fprintf('real: %3.7f \n', t_wr(3,4)); 
    t_w_estimatesr(:, :, i) =t_wr;
    
    t_wo = T_w2e * PB_K * t_c;
    t_wo = [tag_rotation_offset(1:3,1:3) * t_wo(1:3,1:3) t_wo(1:3,4); 0 0 0 1];
    fprintf('opt: %3.7f \n\n', t_wo(3,4)); 
    t_w_estimateso(:, :, i) =t_wo;
    
    zhat_tagr = -1*t_wr(1:3, 3);
    z_tagr =  t_wr(3, 4);
    
    zhat_tago = -1*t_wo(1:3, 3);
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


% load([path, 'dataFromADA/ADA_depthtests_gt.mat'])
% wld_pose = [tags.worldpose(:, :) ones(size(tags.worldpose,1), 1)]';
%%
% 
% init_K
% 
% PB_K
% 
% for i=1:s
%              
%     if(fullcampose)
%         t = tags_train.pose(:,  :, i);
%         t_c = t;       
%     end
%     
%     t_w_init(:, i) =  T_rb2e * init_K * t_c(:,4);
% 
%     t_w_PB(:, i) =  T_rb2e * PB_K * t_c(:,4);
% end
% 
% 
% t_w_init
% 
% t_w_PB
% 
% [mean_abs_error, std_deviation,  mean_error_3D, std_deviation_3D] = translationErrorBetweenPointsInWorld(t_w_init, t_w_PB)
% 

%% Perform World Known Poses Refinement
disp('Perform World Known Poses Refinement')

if(load_depth)
    load([path, 'dataFromADA/ADA_depthtests_gt.mat'])
else
    load([path, 'dataFromADA/ADAtags_14062016_wldPoseTag.mat'])
end

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

tag_in_world_gt = T_w2rb * T_rb2tb * wld_pose ;

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@wldpose_minimization_ADA,x0,[],[],[],[],[],[],[], options, tag_in_world_gt, cam_pose, PB_K, transforms);

% get refined extrinsics
t_est = x(1:3)';

WK_K = [PB_K(1:3, 1:3), t_est ; 0 0 0 1];

%% show results

% matrices
init_K

PB_K

WK_K


% world poses

tag_in_world_gt

t_w_init =  T_rb2e * init_K * cam_pose

t_w_PB =  T_rb2e * PB_K * cam_pose

t_w_WK =  T_rb2e * WK_K * cam_pose  


%% errors 

[mean_abs_error_i, std_deviation_i,  mean_error_3D_i, std_deviation_3D_i] = ...
    translationErrorBetweenPointsInWorld(tag_in_world_gt, t_w_init);

[mean_abs_error_PB, std_deviation_PB,  mean_error_3D_PB, std_deviation_3D_PB] = ...
    translationErrorBetweenPointsInWorld(tag_in_world_gt, t_w_PB);

[mean_abs_error_WK, std_deviation_WK,  mean_error_3D_WK, std_deviation_3D_WK] = ...
    translationErrorBetweenPointsInWorld(tag_in_world_gt, t_w_WK);

t_w_PB2 =  T_rb2e * [init_K(1:3, 1:3) PB_K(1:3, 4); 0 0 0 1] * cam_pose
[mean_abs_error_PB2, std_deviation_PB2,  mean_error_3D_PB2, std_deviation_3D_PB2] = ...
    translationErrorBetweenPointsInWorld(tag_in_world_gt, t_w_PB2)



%% printing

mean_abs_error_i = 1000 * mean_abs_error_i;
fprintf('Init:\n mean : %3.4f , %3.4f , %3.4f \n ', mean_abs_error_i(1), mean_abs_error_i(2), mean_abs_error_i(3));

std_deviation_i = 1000 * std_deviation_i;
fprintf('std  : %3.4f , %3.4f , %3.4f \n ', std_deviation_i(1), std_deviation_i(2), std_deviation_i(3));

mean_error_3D_i = 1000 * mean_error_3D_i;
fprintf('mean : %3.4f \n ', mean_error_3D_i);

std_deviation_3D_i = 1000 * std_deviation_3D_i;
fprintf('std  : %3.4f \n ', std_deviation_3D_i);


mean_abs_error_PB = 1000 * mean_abs_error_PB;
fprintf('\n\n PB  :\n mean : %3.4f , %3.4f , %3.4f \n ', mean_abs_error_PB(1), mean_abs_error_PB(2), mean_abs_error_PB(3));

std_deviation_PB = 1000 * std_deviation_PB;
fprintf('std  : %3.4f , %3.4f , %3.4f \n ', std_deviation_PB(1), std_deviation_PB(2), std_deviation_PB(3));

mean_error_3D_PB = 1000 * mean_error_3D_PB;
fprintf('mean : %3.4f \n ', mean_error_3D_PB);

std_deviation_3D_PB = 1000 * std_deviation_3D_PB;
fprintf('std  : %3.4f \n ', std_deviation_3D_PB);

mean_abs_error_WK = 1000 * mean_abs_error_WK;
fprintf('\n\n WK  :\n mean : %3.4f , %3.4f , %3.4f \n ', mean_abs_error_WK(1), mean_abs_error_WK(2), mean_abs_error_WK(3));

std_deviation_WK = 1000 * std_deviation_WK;
fprintf('std  : %3.4f , %3.4f , %3.4f \n ', std_deviation_WK(1), std_deviation_WK(2), std_deviation_WK(3));

mean_error_3D_WK = 1000 * mean_error_3D_WK;
fprintf('mean : %3.4f \n ', mean_error_3D_WK);

std_deviation_3D_WK = 1000 * std_deviation_3D_WK;
fprintf('std  : %3.4f \n ', std_deviation_3D_WK);




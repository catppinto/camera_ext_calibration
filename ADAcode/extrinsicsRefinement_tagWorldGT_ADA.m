clear
clc
path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/';
addpath([path, 'bkg_code'])
%addpath('/homes/apirespi/cat_workspace/src/ExtrinsicsCalibration/matlab_code/bkg_code')

tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];
%% TRAINING

%% LOADING DATA
load ADA_Data; %% A, K, P
load([path, 'dataFromADA/ADAtags_13062016_wldPoseTag.mat'])

%%  find x

K_wc = T_e2c;

w=rodrigues(K_wc(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[K_wc(1,4), K_wc(2,4), K_wc(3,4)];

wld_pose = [tags.worldpose 1]';
cam_pose = [mean(tags.pose(:, 1:3)) 1]';

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@wldpose_minimization_ADA,x0,[],[],[],[],[],[],[], options, wld_pose, cam_pose, x0);

% get refined extrinsics
t_est = x(1:3)';

KK = [K_wc(1:3, 1:3), t_est ; 0 0 0 1]

K_wc

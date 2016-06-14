
clc
clear 

addpath(genpath('/home/cat/Documents/CMU_Herb/camera_ext_calibration/bkg_code'))
%% choose type of points 

load('DATA/FEUPtags.mat');
load('DATA/FEUP_intrinsics.mat'); 

% c_points = calibration_tags.cam_points;
% 
uv_points = (calibration_tags.corners.c1 + calibration_tags.corners.c2 + calibration_tags.corners.c3 + calibration_tags.corners.c4)/4
uv_points = [uv_points, ones(size(uv_points,1),1)];       

wld_points = calibration_tags.wld_points;
wld_points = [wld_points, ones(size(wld_points,1),1)];
wld_points(1:4, 3) = 0.005*rand(4,1)

%% calculte extrinsics : first step

c = CamExtrinsicsCalibrationClass();

c = c.CalculateExtrinsics_wldcampointsIntrinsics(wld_points, uv_points, intrinsics)
            

%% function evaluation 

[res_lsq_noJacobian, res_lsq_Jacobian] = camCalib_fixedIntrinsics_lsqcurvefit_func(c)

%%  displays
A = c.intrinsic_camera_matrix;
disp(' Ground Truth ') 
c.extrinsic_camera_matrix

disp(' Initial Guess ') 
K_svd = [c.R_ext c.t_ext; 0 0 0 1]
    
disp(' LM no Jacobian ') 
K_lm_noJac = res_lsq_noJacobian.ext_est

disp(' LM Jacobian on ') 
K_lm_Jac = res_lsq_Jacobian.ext_est

%% TESTING

uv_points

uv_est_svd = (intrinsics * K_svd * wld_points')';
% (inv(intrinsics * K_svd) * uv_est_svd')'
% (inv(intrinsics * K_svd) * uv_points')'
% wld_points
uv_est_svd = [uv_est_svd(:, 1)./uv_est_svd(:, 3) uv_est_svd(:, 2)./uv_est_svd(:, 3) uv_est_svd(:, 3)./uv_est_svd(:, 3)]

uv_est_lm_noJac = (intrinsics * K_lm_noJac * wld_points')';
uv_est_lm_noJac = [uv_est_lm_noJac(:, 1)./uv_est_lm_noJac(:, 3) uv_est_lm_noJac(:, 2)./uv_est_lm_noJac(:, 3) uv_est_lm_noJac(:, 3)./uv_est_lm_noJac(:, 3)]

uv_est_lm_Jac = (intrinsics * K_lm_Jac * wld_points')';
uv_est_lm_Jac = [uv_est_lm_Jac(:, 1)./uv_est_lm_Jac(:, 3) uv_est_lm_Jac(:, 2)./uv_est_lm_Jac(:, 3) uv_est_lm_Jac(:, 3)./uv_est_lm_Jac(:, 3)]

%% Reprojection Error 




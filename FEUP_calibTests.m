clc
clear 

addpath(genpath('/home/apirespi/Documents/Thesis/ExtrinsicsCalibration/matlab_code/bkg_code'))
addpath('/home/cat/Documents/CMU_Herb/camera_ext_calibration/FEUPcalibration')
%% choose type of points 
use_random_points = false;

[cam_points wld_points]= getDatafromFEUP();

% n = size(wld_points, 1);
% max_z = 0.93;
% min_z = 0.90;
% wz = min_z + (max_z-min_z).*rand(n,1);
% 
% wld_points(:, 3) = wz;

        
intrinsics = [ 529.6085148771313, 0, 318.3732583194404; ...
            0, 529.3631821305655, 243.8970374362488; ...
            0, 0, 1];
            


%% calculte extrinsics : first step

c = CamExtrinsicsCalibrationClass();

c = c.CalculateExtrinsics_wldcampointsIntrinsics(wld_points, cam_points, intrinsics);

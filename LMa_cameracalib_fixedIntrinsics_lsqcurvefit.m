
clc
clear 

addpath(genpath('/home/apirespi/Documents/Thesis/ExtrinsicsCalibration/matlab_code/bkg_code'))

%% choose type of points 
use_random_points = false;

if(use_random_points)
    n = 100;
    wx = 0 + (2-0).*rand(n,1); 
    wy = 0 + (2-0).*rand(n,1); 
    wz = 0 + (5-0).*rand(n,1);

    wld_points = [wx wy wz ones(n, 1)];
else
    load test1_data.mat
    table_file1 = tags_file1(1);
    wld_points = table_file1.wld(1:3);
    cam_points = table_file1.cam_middlep;
    
end

%% calculte extrinsics : first step

c = CamExtrinsicsCalibrationClass();

if(use_random_points)
    c = c.CalculateExtrinsics(wld_points, true);
else
    c = c.CalculateExtrinsics_wldcampoints(wld_points, cam_points);
end

if(use_random_points)
    %% function evaluation 

    [res_lsq_noJacobian, res_lsq_Jacobian] = camCalib_fixedIntrinsics_lsqcurvefit_func(c)

    %%  displays
    K = c.intrinsic_camera_matrix;
    disp(' Ground Truth ') 
    c.extrinsic_camera_matrix

    disp(' Initial Guess ') 
    m = res_lsq_noJacobian.ext0

    disp(' LM no Jacobian ') 
    m = res_lsq_noJacobian.ext_est
    p = K * m;
    error_1 = computeReprojectionError( c.wld_points, c.proj_matrix, p);

    disp(' LM Jacobian on ') 
    m = res_lsq_Jacobian.ext_est
    p = K * m;
    error_2 = computeReprojectionError( c.wld_points, c.proj_matrix, p);

    fprintf('Error no JACOBIAN : %f \n Error w JACOBIAN : %f \n', error_1, error_2);  
end

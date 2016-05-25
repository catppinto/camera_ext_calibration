%%% TODO : verificar a reprojection dos pontos e como estar a ser calculado
%%% o erro. Estou a usar os valores da camera de de u,v,w qd so devia ser
%%% preciso u,v. Isto esta a ocorrer na funcao  computeReprojectionError.
%%% Deve-se dividir por w? Verificar isto ... 
%%% TODO : testar com rand para o barulho

addpath('bkg_code')

clc
clear 


c = CamExtrinsicsCalibrationClass();

n = 100;
wx = 0 + (2-0).*rand(n,1); 
wy = 0 + (2-0).*rand(n,1); 

min_z = 0;
max_z = 2;
wz = min_z + (max_z-min_z).*rand(n,1);

wld_points = [wx wy wz ones(n, 1)];

%%
c = c.CalculateExtrinsics(wld_points, true);

%% Error calculation 


disp('Extrinsics GT');
c.extrinsic_camera_matrix

disp('Extrinsics Comp');
ext_comp = [c.R_ext, c.t_ext; 0 0 0 1]

proj_matrix_fromExtrinsicCalculus = c.intrinsic_camera_matrix * ext_comp;
reprojection_error = computeReprojectionError( c.wld_points, c.proj_matrix, proj_matrix_fromExtrinsicCalculus) 

xyzCam_gt = (c.extrinsic_camera_matrix * c.wld_points')' ; 
xyzCam_test = (ext_comp * c.wld_points')' ; 
error_xyzCam = sum(sum((xyzCam_gt - xyzCam_test).^2))


%% 

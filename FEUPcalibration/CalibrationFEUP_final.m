% close all
% clc
% clear all

n= 4;

w = [0.7005 0.024 0; ... 
   0.4 0.024 0 ; ...
   0.2495 0.024 0; ...
   0.402 0.19 0; ...
   0.449 -0.004 0.16627; ...
   0.449 -0.004 0.11877; ...
   0.449 -0.004 0.07127; ...
   0.449 -0.004 0.02377];

ctag = [-0.322713 0.0119508 0.671776; ...
    0.0577349 0.0274597 0.670989; ...
    0.243119 0.0353418 0.669769; ...
    0.053183 0.158643 0.511146; ...
    0.0110721 -0.161089 0.57441; ...
    0.00809416 -0.114525 0.611413; ...
    0.00515552 -0.0678623 0.646763; ...
    0.00249965 -0.0216926 0.685821];
ctag_h = [ctag ones(8,1)]';

K_rigid = performCalibration(w, ctag); 

%% point in world frame
[w_p] = getWorldPointFromCameraPoint(ctag, K_rigid);
w_test = w_p(n, :); 

%%  find x
  
x0 =[K_rigid(1,4), K_rigid(2,4), K_rigid(3,4)];

ctag = [ctag ones(size(ctag,1), 1)]'; 
w = [w ones(size(w,1), 1)]';
%% World Known Poses Based Optimization
options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@wldpose_minimization_FEUP,x0,[],[],[],[],[],[],[], options, w, ctag, K_rigid);

% get refined extrinsics
t_est = x(1:3)';

K_WK = [K_rigid(1:3, 1:3), t_est ; 0 0 0 1];

%% Results

% matrices

K_rigid

K_WK

% world poses

t_w_rigid =  K_rigid * ctag

t_w_wk =  K_WK * ctag 


%% error of projection 

[mean_abs_error_rig, std_deviation_rig,  mean_error_3D_rig, std_deviation_3D_rig] = translationErrorBetweenPointsInWorld(w, t_w_rigid)

[mean_abs_error_wk, std_deviation_wk,  mean_error_3D_wk, std_deviation_3D_wk]  = translationErrorBetweenPointsInWorld(w, t_w_wk)







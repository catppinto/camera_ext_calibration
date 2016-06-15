clear
clc
load('DATA/FEUPtags.mat');
load('DATA/FEUPdata.mat'); 
addpath('/home/cat/Documents/CMU_Herb/camera_ext_calibration/bkg_code')
addpath('/home/cat/Documents/CMU_Herb/camera_ext_calibration/ADAcode')
%% LOADING DATA

%%  find x

K_ec = [   -0.9988   -0.0278    0.0403 0.4249;
   -0.0488    0.6135   -0.7882 0.5380;
   -0.0028   -0.7892   -0.6141 0.4216;
   0 0 0 1];
  
w=rodrigues(K_ec(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[K_ec(1,4), K_ec(2,4), K_ec(3,4)];

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
ctag = [ctag ones(size(ctag,1), 1)]'; 
w = [w ones(size(w,1), 1)]';
%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@wldpose_minimization_FEUP,x0,[],[],[],[],[],[],[], options, w, ctag, K_ec);

% get refined extrinsics
t_est = x(1:3)';

k_ec_opt = [K_ec(1:3, 1:3), t_est ; 0 0 0 1]

K_ec

t_w_opt =  k_ec_opt * ctag 

t_w_prior_opt =  K_ec * ctag

w


%% error of projection 

[error_t_opt, overall_opt] = translationErrorBetweenPointsInWorld(w, t_w_opt)
[error_t_prioropt, overall_prior_opt ]= translationErrorBetweenPointsInWorld(w, t_w_prior_opt)







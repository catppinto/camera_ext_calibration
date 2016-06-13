%%%% TODO 
% add Intrinsics from ADA
% add corners (redo calibration)
% Confirm computation of Tw
%%%%

clear
clc
path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/';
addpath([path, 'bkg_code'])
%addpath('/homes/apirespi/cat_workspace/src/ExtrinsicsCalibration/matlab_code/bkg_code')

tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];
%% TRAINING

%% LOADING DATA
load ADA_Data; %% A, K, P
load([path, 'dataFromADA/ADAtags_13062016_trueSizeTag.mat'])
corners = [];

%%  find x

K_wc = T_e2c;

w=rodrigues(K_wc(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4)]

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@minimize_tableTagErrorADA,x0,[],[],[],[],[],[],[], options, tags_train.pose, tags_train.corners, x0);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

K_wc

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
for i=1:size(tags_train.pose,1) 
     
    t = tags_train.pose(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    t_wr = T_tb2rb * T_rb2e * K_wc * t_c;
    t_wr = tag_rotation_offset * t_wr;
    fprintf('real: %3.7f \n', t_wr(3,4)); 
    t_w_estimatesr(:, :, i) =t_wr;
    
    t_wo = T_tb2rb * T_rb2e * KK * t_c;
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

% %% corners
% t_w_estimates = t_w_estimatesr;
% tag_size = 0.06;
% ts = tag_size/2;
% for i=1:size(t_w_estimates,3)
%     t = t_w_estimates(:, :, i);
%     corner1_est = t; corner1_est(1,4) = t(1,4)+ts; corner1_est(2,4) = t(2,4)+ts; 
%     corner2_est = t; corner2_est(1,4) = t(1,4)-ts; corner2_est(2,4) = t(2,4)+ts;
%     corner3_est = t; corner3_est(1,4) = t(1,4)-ts; corner3_est(2,4) = t(2,4)-ts;
%     corner4_est = t; corner4_est(1,4) = t(1,4)+ts; corner4_est(2,4) = t(2,4)-ts;
%     
%     T_w1(:, :, i) = corner1_est; 
%     T_w2(:, :, i) = corner2_est;
%     T_w3(:, :, i) = corner3_est;
%     T_w4(:, :, i) = corner4_est;
% end
% 
% for n =1:size(t_w_estimates,3)
%    T_c1(:, :, n) = inv(K_wc) * T_w1(:, :, n);
%    T_c2(:, :, n) = inv(K_wc) * T_w2(:, :, n);
%    T_c3(:, :, n) = inv(K_wc) * T_w3(:, :, n);
%    T_c4(:, :, n) = inv(K_wc) * T_w4(:, :, n);
% end
%         
% c1_est = []; c2_est= []; c3_est= []; c4_est= [];
% for n =1:size(t_w_estimates,3)
%    tc_1 = T_c1(:, :, n) ;
%    corner1_est = A * tc_1(1:4, 4);
%    c1_est(n, :) = [corner1_est(1)/corner1_est(3) corner1_est(2)/corner1_est(3) 1];
% 
%    tc_2 = T_c2(:, :, n) ;
%    corner2_est = A * tc_2(1:4, 4);
%    c2_est(n, :) = [corner2_est(1)/corner2_est(3) corner2_est(2)/corner2_est(3) 1];
%    
%    tc_3 = T_c3(:, :, n) ;
%    corner3_est = A * tc_3(1:4, 4);
%    c3_est(n, :) = [corner3_est(1)/corner3_est(3) corner3_est(2)/corner3_est(3) 1];
%    
%    tc_4 = T_c4(:, :, n) ;
%    corner4_est = A * tc_4(1:4, 4);
%    c4_est(n, :) = [corner4_est(1)/corner4_est(3) corner4_est(2)/corner4_est(3) 1];
% end
% 
% for n =1:size(t_w_estimates,3)
%     norm_c1(n) = norm(c1_est(n, :) - corners.c3(n, :));
%     norm_c2(n) = norm(c2_est(n, :) - corners.c4(n, :));
%     norm_c3(n) = norm(c3_est(n, :) - corners.c1(n, :));
%     norm_c4(n) = norm(c4_est(n, :) - corners.c2(n, :));
% end
% 
% reprojection_error = sum(norm_c1) + sum(norm_c2) + sum(norm_c3) + sum(norm_c4);
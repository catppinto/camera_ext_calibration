addpath(strcat(pwd, '/bkg_code'))

clear
clc
%% load data
load test1_data.mat;
% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P
load HERB_tableData;

%% LOADING DATA

marker_wc_q = tags_file1(1).wld;
marker_wc_rot = quatToRotationMatrix(marker_wc_q(4:7));
marker_wc_t = marker_wc_q(1:3)';
Xc = [marker_wc_rot marker_wc_t ; 0 0 0 1]

marker_rot_tabl = (tb_kinbody_offset(1:3, 1:3));
%Xc_ = [inv(marker_rot_tabl)*marker_wc_rot marker_wc_t ; 0 0 0 1]
%Xw = K_wc * Xc_ 
Xw = K_wc * Xc

% uv point
%uv = A * Xc(1:4, 4);
uv = tags_file1(1).cam_middlep;
uv = [uv(1)/uv(3) uv(2)/uv(3) 1]'

%% find x0
tag_x = Xw(1,4)
tag_y = Xw(2,4) 
r_ = atan2(Xw(1,1), Xw(2,1))

w=rodrigues(K_wc(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4), r_, tag_x, tag_y]

%% NOT NEEDED : JUST TO SEE THE DATA
x = x0;
tb_height = 0.7380;
observed_proj_marker = uv;
                        
R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

KK_wc = [R_wc t_wc ; 0 0 0 1];

r_ = x(7);
tb_pose = x(8:9);
expected_marker = [cos(r_) -sin(r_)  0   tb_pose(1) ; ...
                    sin(r_) cos(r_)  0   tb_pose(2) ; ...
                    0       0        1   tb_height; ...
                    0       0        0   1];

                
marker_c = inv(KK_wc) * expected_marker;
expected_proj_marker_Xc = marker_c(:, 4);

expected_proj_marker_uv = A * expected_proj_marker_Xc;

observed_proj_marker_uv = [observed_proj_marker(1)./observed_proj_marker(3); ...
                        observed_proj_marker(2)./observed_proj_marker(3); ...
                        1]
expected_proj_marker_uv = [expected_proj_marker_uv(1)./expected_proj_marker_uv(3); ...
                        expected_proj_marker_uv(2)./expected_proj_marker_uv(3); ...
                        1]
                           
F = norm(observed_proj_marker_uv - expected_proj_marker_uv);                    
                    
% NOT NEEDED : UNTIL HERE
%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');
lb = [-Inf, -Inf, -Inf, -5, -5, -5, -Inf, -5, -5];
ub = [Inf, Inf, Inf, 5, 5, 5, Inf, 5, 5];
[x, fval, exittag] = fmincon(@find_expectedonetagpoint_projection,x0,[],[],[],[],lb,ub,[], options);


% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

KK_wc

% get refined tb pose
rxy_table = x(7:9) ;
r = rxy_table(1); x=rxy_table(2); y=rxy_table(3);
tb_pose_result = [ cos(r), -sin(r), 0 , x; ...
            sin(r), cos(r), 0, y ; ...
            0, 0, 1, tb_height; ...
            0 0 0 1]
             
tb_pose = [ cos(r_), -sin(r_), 0 , Xw(1); ...
    sin(r_), cos(r_), 0, Xw(2); ...
    0, 0, 1, Xw(3); ...
    0 0 0 1]

error_tb_position = sqrt( (tb_pose(1, 4)-tb_pose_result(1, 4))^2 + (tb_pose(2, 4)-tb_pose_result(2, 4))^2 + (tb_pose(3, 4)-tb_pose_result(3, 4))^2)
error_t=  sqrt( (KK(1, 4)-KK_wc(1, 4))^2 + (KK(2, 4)-KK_wc(2, 4))^2 + (KK(3, 4)-KK_wc(3, 4))^2)

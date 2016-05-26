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
Xc_ = [inv(marker_rot_tabl)*marker_wc_rot marker_wc_t ; 0 0 0 1]
Xw = K * Xc_ 
uv = P * Xw
uv = [uv(1,4)/uv(3,4) uv(2,4)/uv(3,4) 1]'

%% find x0
tag_x = Xw(1,4)
tag_y = Xw(2,4) 
r_ = atan2(Xw(1,1), Xw(2,1))

w=rodrigues(K(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K(1,4), K(2,4), K(3,4), r_, tag_x, tag_y]

%% NOT NEEDED : JUST TO SEE THE DATA
x = x0;
tb_height = 0.7380;
observed_proj_marker = uv;
     
R = rodrigues(x(1:3));
t=[x(4);x(5);x(6)];

KK = [R t ; 0 0 0 1];
PP = A * KK ;

r_ = x(7);
tb_pose = x(8:9);
expected_marker = [cos(r_) -sin(r_)  0   tb_pose(1) ; ...
                    sin(r_) cos(r_)  0   tb_pose(2) ; ...
                    0       0        1   tb_height; ...
                    0       0        0   1]
marker_cam = PP * expected_marker;
expected_proj_marker = marker_cam(1:3, 4);
observed_proj_marker = [observed_proj_marker(1)./observed_proj_marker(3); ...
                        observed_proj_marker(2)./observed_proj_marker(3); ...
                        1]
expected_proj_marker = [expected_proj_marker(1)./expected_proj_marker(3); ...
                        expected_proj_marker(2)./expected_proj_marker(3); ...
                        1]
                    
% NOT NEEDED : UNTIL HERE
%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');
lb = [-Inf, -Inf, -Inf, -5, -5, -5, -Inf, -5, -5];
ub = [Inf, Inf, Inf, 5, 5, 5, Inf, 5, 5];
[x, fval, exittag] = fmincon(@find_expectedmarker_projection,x0,[],[],[],[],lb,ub,[], options);


% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

K

% get refined tb pose
rxy_table = x(7:9) ;
r = rxy_table(1); x=rxy_table(2); y=rxy_table(3);
tb_pose_result = [ cos(r), -sin(r), 0 , x; ...
            sin(r), cos(r), 0, y ; ...
            0, 0, 1, tb_height; ...
            0 0 0 1]
             
tb_pose = [ cos(r_), -sin(r_), 0 , Xw(1,4); ...
    sin(r_), cos(r_), 0, Xw(2,4); ...
    0, 0, 1, Xw(3,4); ...
    0 0 0 1]

error_tb_position = sqrt( (tb_pose(1, 4)-tb_pose_result(1, 4))^2 + (tb_pose(2, 4)-tb_pose_result(2, 4))^2 + (tb_pose(3, 4)-tb_pose_result(3, 4))^2)
error_t=  sqrt( (KK(1, 4)-K(1, 4))^2 + (KK(2, 4)-K(2, 4))^2 + (KK(3, 4)-K(3, 4))^2)

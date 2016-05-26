clear
clc
% load data
load test1_data.mat;
load HERB_projectionData; %% A, K, P
load HERB_tableData;

% LOADING DATA

marker_wc_q = tags_file1(1).wld;
marker_wc_rot = quatToRotationMatrix(marker_wc_q(4:7));
marker_wc_t = marker_wc_q(1:3)';
Xc = [marker_wc_rot marker_wc_t ; 0 0 0 1]

marker_rot_tabl = (tb_kinbody_offset(1:3, 1:3));
Xc_ = [inv(marker_rot_tabl)*marker_wc_rot marker_wc_t ; 0 0 0 1]
Xw = K * Xc_ 
uv = P * Xw
uv = [uv(1,4)/uv(3,4) uv(2,4)/uv(3,4) 1]'

tb_height = 0.7380;
    
tag_x = Xw(1,4);
tag_y = Xw(2,4);
r_ = atan2(Xw(1,1), Xw(2,1));

w=rodrigues(K(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x =[wx, wy, wz, K(1,4), K(2,4), K(3,4), r_, tag_x, tag_y];

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
expected_proj_marker = [expected_proj_marker(1)./expected_proj_marker(3); ...
                        expected_proj_marker(2)./expected_proj_marker(3); ...
                        1]
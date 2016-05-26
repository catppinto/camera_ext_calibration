function [F] = find_expectedmarker_projection(x)

%% load data
load test1_data.mat;

%% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P

load HERB_tableData;

%% get marker data 

marker_wc_q = tags_file1(1).wld;
marker_wc_rot = quatToRotationMatrix(marker_wc_q(4:7));
marker_wc_t = marker_wc_q(1:3)';
Xc = [marker_wc_rot marker_wc_t ; 0 0 0 1];

marker_rot_tabl = (tb_kinbody_offset(1:3, 1:3));
Xc_ = [inv(marker_rot_tabl)*marker_wc_rot marker_wc_t ; 0 0 0 1];
Xw = K * Xc_ ;
uv = P * Xw;
uv = [uv(1,4)/uv(3,4) uv(2,4)/uv(3,4) 1]';

%% Observed marker

observed_proj_marker = uv; 
     

%% Expected marker

R = rodrigues(x(1:3));
t=[x(4);x(5);x(6)];

KK = [R t ; 0 0 0 1];
PP = A * KK ;

r_ = x(7);
tb_pose = x(8:9);
expected_marker = [cos(r_) -sin(r_)  0   tb_pose(1) ; ...
                    sin(r_) cos(r_)  0   tb_pose(2) ; ...
                    0       0        1   tb_height; ...
                    0       0        0   1];

                
marker_cam = PP * expected_marker;
expected_proj_marker = marker_cam(1:3, 4);

observed_proj_marker = [observed_proj_marker(1)./observed_proj_marker(3); ...
                        observed_proj_marker(2)./observed_proj_marker(3); ...
                        1];
expected_proj_marker = [expected_proj_marker(1)./expected_proj_marker(3); ...
                        expected_proj_marker(2)./expected_proj_marker(3); ...
                        1];
                           
F = norm(observed_proj_marker - expected_proj_marker);





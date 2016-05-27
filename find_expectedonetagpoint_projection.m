function [F] = find_expectedonetagpoint_projection(x)

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
Xw = K_wc * Xc_ ;
uv = P * Xw(:, 4);
uv = [uv(1)/uv(3) uv(2)/uv(3) 1]';

%% Observed marker

observed_proj_marker = uv; 
     

%% Expected marker

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
                        1];
expected_proj_marker_uv = [expected_proj_marker_uv(1)./expected_proj_marker_uv(3); ...
                        expected_proj_marker_uv(2)./expected_proj_marker_uv(3); ...
                        1];
                           
F = norm(observed_proj_marker_uv - expected_proj_marker_uv);





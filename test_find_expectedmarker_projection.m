function [F] = test_find_expectedmarker_projection(x)


%% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P
load 'test_data_extrinsicsTagOnly.mat';

%% 

tb_height = 0.7380;


%% Observed marker
marker_u = cam_points(1:3,4);
observed_proj_marker = marker_u;
     
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

%%
observed_proj_marker = [observed_proj_marker(1)./observed_proj_marker(3); ...
                        observed_proj_marker(2)./observed_proj_marker(3); ...
                        1];
expected_proj_marker = [expected_proj_marker(1)./expected_proj_marker(3); ...
                        expected_proj_marker(2)./expected_proj_marker(3); ...
                        1];

                           
F = norm(observed_proj_marker - expected_proj_marker);





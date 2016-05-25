function [F] = find_expectedmarker_projection(x)

%% load data
load test1_data.mat;

%% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P

load HERB_tableData;

%% get marker data 

% better ideia to use marker from cam ??? 
%marker_cam = tags_file1(1).cam_middlep; 

marker_uv = tags_file1(1).cam_middlep;
marker_wc_q = tags_file1(1).wld;
marker_wc_rot = quatToRotationMatrix(marker_wc_q(4:7));
marker_wc_t = marker_wc_q(1:3)';
marker_wc = [marker_wc_rot marker_wc_t ; 0 0 0 1];

%% Observed marker

%observed marker
%observed_proj_marker = P * marker_w;
%observed_proj_marker = observed_proj_marker(1:3)';

% using camera points
observed_proj_marker = marker_uv'; 
     


%% Expected marker

% convert the 3-vector [wx wy wz] of the Rodrigues representation
% into the 3x3 rotation matrix
% theta2=x(1)^2+x(2)^2+x(3)^2
% theta=sqrt(theta2)
% omega=[0 -x(3) x(2);
%     x(3) 0  -x(1);
%     -x(2) x(1) 0;];
% R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/(theta^2))*(omega*omega)

R = rodrigues(x(1:3));
t=[x(4);x(5);x(6)];

K = [R t ; 0 0 0 1];
P = A * K ;

r_ = x(7);
tb_pose = x(8:9);
expected_marker = [cos(r_) -sin(r_)  0   tb_pose(1) ; ...
                    sin(r_) cos(r_)  0   tb_pose(2) ; ...
                    0       0        1   tb_height; ...
                    0       0        0   1];
marker_cam = P * expected_marker;
expected_proj_marker = marker_cam(1:3, 4);
expected_proj_marker = [expected_proj_marker(1)./expected_proj_marker(3); ...
                        expected_proj_marker(2)./expected_proj_marker(3); ...
                        1];

%gilwoo version 
% expected_proj_marker_gilwoo_ = inv(K) * expected_marker * inv(tb_kinbody_offset);
% expected_proj_marker_gilwoo = P * expected_proj_marker_gilwoo_(:,4)
% 
% expected_proj_marker_gilwoo = [expected_proj_marker_gilwoo(1)./expected_proj_marker_gilwoo(3); ...
%                                expected_proj_marker_gilwoo(2)./expected_proj_marker_gilwoo(3); ...
%                                1]

                           
F = norm(observed_proj_marker - expected_proj_marker);





function [F] = minimize_extrinsicAndPose(x)

%% LOADING DATA
[tagcorners, halftagsize] = loadTagCornersdata('test1_data.mat') ;

%% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P

load HERB_tableData;

%% Observed marker

observed_uv_c1 = tagcorners.uv.c1'; 
observed_uv_c2 = tagcorners.uv.c2'; 
observed_uv_c3 = tagcorners.uv.c3'; 
observed_uv_c4 = tagcorners.uv.c4'; 

%% Conversion to R,t extrinsiccs and Tag pose 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

KK_wc = [R_wc t_wc ; 0 0 0 1];

tag_theta = x(7);
tag_x = x(8);
tag_y = x(9);

%% Expected world coordinates calculus 

expected_rotation = [   cos(tag_theta) -sin(tag_theta)  0 ; ...
                        sin(tag_theta)  cos(tag_theta)  0 ; ...
                        0                            0  1];
expected_w_c1 = [tag_x+halftagsize, tag_y-halftagsize, tb_height]';
expected_w_c2 = [tag_x+halftagsize, tag_y+halftagsize, tb_height]';
expected_w_c3 = [tag_x-halftagsize, tag_y+halftagsize, tb_height]';
expected_w_c4 = [tag_x-halftagsize, tag_y-halftagsize, tb_height]';

expected_c1 = [ expected_rotation expected_w_c1; 0 0 0 1];
expected_c2 = [ expected_rotation expected_w_c2; 0 0 0 1];
expected_c3 = [ expected_rotation expected_w_c3; 0 0 0 1];
expected_c4 = [ expected_rotation expected_w_c4; 0 0 0 1];
           
% corner 1 
expected_c1_c = inv(KK_wc) * expected_c1;
expected_c1_Xc = expected_c1_c(:, 4);
expected_uv_c1 = A * expected_c1_Xc;
expected_uv_c1 = [expected_uv_c1(1)/expected_uv_c1(3), expected_uv_c1(2)/expected_uv_c1(3) 1]';

% corner 2
expected_c2_c = inv(KK_wc) * expected_c2;
expected_c2_Xc = expected_c2_c(:, 4);
expected_uv_c2 = A * expected_c2_Xc;
expected_uv_c2 = [expected_uv_c2(1)/expected_uv_c2(3), expected_uv_c2(2)/expected_uv_c2(3) 1]';

% corner 3 
expected_c3_c = inv(KK_wc) * expected_c3;
expected_c3_Xc = expected_c3_c(:, 4);
expected_uv_c3 = A * expected_c3_Xc;
expected_uv_c3 = [expected_uv_c3(1)/expected_uv_c3(3), expected_uv_c3(2)/expected_uv_c3(3) 1]';

% corner 4 
expected_c4_c = inv(KK_wc) * expected_c4;
expected_c4_Xc = expected_c4_c(:, 4);
expected_uv_c4 = A * expected_c4_Xc;
expected_uv_c4 = [expected_uv_c4(1)/expected_uv_c4(3), expected_uv_c4(2)/expected_uv_c4(3) 1]';

%% cost function 
norm_c1 = norm(observed_uv_c1 - expected_uv_c1);
norm_c2 = norm(observed_uv_c2 - expected_uv_c2);
norm_c3 = norm(observed_uv_c3 - expected_uv_c3);
norm_c4 = norm(observed_uv_c4 - expected_uv_c4);

F = norm_c1 + norm_c2 + norm_c3 + norm_c4 ;





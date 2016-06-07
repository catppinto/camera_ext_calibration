addpath(strcat(pwd, '/bkg_code'))

clear
clc
%% load data
load HERB_projectionData; %% A, K, P
load HERB_tableData;

%% LOADING DATA

[tagcorners, halftagsize] = loadTagCornersdata('test1_data.mat') ;

%%  find x0

tag_x = tagcorners.xyz_w(1);
tag_y = tagcorners.xyz_w(2);
tag_theta = atan2(tagcorners.rot_w(1,1), tagcorners.rot_w(2,1));
tag_ = [tagcorners.rot_w tagcorners.xyz_w; 0 0 0 1];

w=rodrigues(K_wc(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4), tag_theta, tag_x, tag_y]

 %% NOT NEEDED : JUST TO SEE THE DATA
x = x0;
observed_uv_c1 = tagcorners.uv.c1'; 
observed_uv_c2 = tagcorners.uv.c2'; 
observed_uv_c3 = tagcorners.uv.c3'; 
observed_uv_c4 = tagcorners.uv.c4'; 

% Conversion to R,t extrinsiccs and Tag pose 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

KK_wc = [R_wc t_wc ; 0 0 0 1];

tag_theta = x(7);
tag_x = x(8);
tag_y = x(9);

% Expected world coordinates calculus 

expected_rotation = [   cos(tag_theta) -sin(tag_theta)  0 ; ...
                        sin(tag_theta)  cos(tag_theta)  0 ; ...
                        0                            0  1]
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

observed_uv_c1
expected_uv_c1

observed_uv_c2
expected_uv_c2

observed_uv_c3
expected_uv_c3

observed_uv_c4
expected_uv_c4


% cost function 
norm_c1 = norm(observed_uv_c1 - expected_uv_c1);
norm_c2 = norm(observed_uv_c2 - expected_uv_c2);
norm_c3 = norm(observed_uv_c3 - expected_uv_c3);
norm_c4 = norm(observed_uv_c4 - expected_uv_c4);                
                     
% NOT NEEDED : UNTIL HERE
%% optimization using fmincon


options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');
lb = [-Inf, -Inf, -Inf, -5, -5, -5, -Inf, -5, -5];
ub = [Inf, Inf, Inf, 5, 5, 5, Inf, 5, 5];
[x, fval, exittag] = fmincon(@minimize_extrinsicAndPose,x0,[],[],[],[],lb,ub,[], options);


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
             
tb_pose = tag_

error_tb_position = sqrt( (tb_pose(1, 4)-tb_pose_result(1, 4))^2 + (tb_pose(2, 4)-tb_pose_result(2, 4))^2 + (tb_pose(3, 4)-tb_pose_result(3, 4))^2)
error_t=  sqrt( (KK(1, 4)-KK_wc(1, 4))^2 + (KK(2, 4)-KK_wc(2, 4))^2 + (KK(3, 4)-KK_wc(3, 4))^2)

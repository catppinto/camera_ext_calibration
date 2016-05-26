clear
clc

addpath(strcat(pwd, '/bkg_code'))
load 'HERB_projectionData.mat' % P, A, K

%% load data

n = 1;
min_x = 0.5; max_x = 0.7;
min_y = -0.3; max_y = 0.3;
min_z = 0.74; max_z = 0.76;
wx = min_x + (max_x-min_x).*rand(n,1); 
wy = min_y + (max_y-min_y).*rand(n,1); 
wz = min_z + (max_z-min_z).*rand(n,1);
wld_points = [wx wy wz ones(n, 1)];

a = 0.8 + (1-0.8).*rand(n,1); 
e = 0.8 + (1-0.8).*rand(n,1); 
i = 0.8 + (1-0.8).*rand(n,1); 
b = 0 + (0.3-0).*rand(n,1); 
c = 0 + (0.3-0).*rand(n,1); 
d = 0 + (0.3-0).*rand(n,1); 
f = 0 + (0.3-0).*rand(n,1); 
g = 0 + (0.3-0).*rand(n,1); 
h = 0 + (0.3-0).*rand(n,1); % not the best approximation in the world since is not a rotation matrix per si

rot_matrix = [a b c; d e f; g h i];
p = [ [rot_matrix ;[0 0 0]], wld_points']

cam_points = P*p;
cam_points_prior_noise = cam_points;

%add noise
snr = 50;
cam_points = [awgn(cam_points(1:3, :), snr); 0 0 0 1];

save('test_data_extrinsicsTagOnly.mat', 'cam_points', 'wld_points');

%% LOADING DATA

marker_w = wld_points(1:3)';
marker_u = cam_points(1:3,4);

%% find x0
tag_x = marker_w(1);
tag_y = marker_w(2);
r_ = atan2(p(1,1),p(2,1)) ;

w=rodrigues(K(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K(1,4), K(2,4), K(3,4), r_, tag_x, tag_y];

%% NOT NEEDED : JUST TO SEE THE DATA
x = x0;
tb_height = 0.7380;
marker_u = cam_points(1:3,4);
observed_proj_marker = marker_u;
     
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
                    
%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');
lb = [-Inf, -Inf, -Inf, -5, -5, -5, -Inf, -5, -5];
ub = [Inf, Inf, Inf, 5, 5, 5, Inf, 5, 5];
[x, fval, exittag] = fmincon(@test_find_expectedmarker_projection,x0,[],[],[],[],lb,ub,[], options);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

K

% get refined tb pose
rxy_table = x(7:9) ;
r = rxy_table(1); x=rxy_table(2); y=rxy_table(3);
tb_pose_result = [ cos(r), -sin(r), 0 , x; ...
            sin(r), cos(r), 0, y ; ...
            0, 0, 1, 0.738; ...
            0 0 0 1]
             
tb_pose = [ cos(r_), -sin(r_), 0 , marker_w(1); ...
    sin(r_), cos(r_), 0, marker_w(2); ...
    0, 0, 1, marker_w(3); ...
    0 0 0 1]
        
error_tb_position = sqrt( (tb_pose(1, 4)-tb_pose_result(1, 4))^2 + (tb_pose(2, 4)-tb_pose_result(2, 4))^2 + (tb_pose(3, 4)-tb_pose_result(3, 4))^2)
error_t=  sqrt( (KK(1, 4)-K(1, 4))^2 + (KK(2, 4)-K(2, 4))^2 + (KK(3, 4)-K(3, 4))^2)
       
        
        



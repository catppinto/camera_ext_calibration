addpath(strcat(pwd, '/bkg_code'))

%% load data
load test1_data.mat;
% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P

load HERB_tableData;

%% LOADING DATA

marker_wc_q = tags_file1(1).wld;
marker_wc_rot = quatToRotationMatrix(marker_wc_q(4:7));
marker_wc_rot = [ 0.3027         0    0.0004; -0.0002    0.2772    0.4752;  -0.0001   -0.4752    0.2772];

marker_wc_t = marker_wc_q(1:3)';
marker_wc = [marker_wc_rot marker_wc_t ; 0 0 0 1];

%% obs marker + exp marker
% obs marker : in optimization function
% expected marker : 

% gilwoo : 
%table_pose = K * marker_w * tb_kinbody_offset;
%r_ = atan2(table_pose(1,1), table_pose(2,1));

% convert from xyz from camera to xyz from world, using first estimates of
% R and t external

marker_rot_tabl = [tb_kinbody_offset(1:3, 1:3) [0;0;0] ; 0 0 0 1];
marker_wc = inv(marker_rot_tabl) * marker_wc
marker_w = K * marker_wc % TODO confirm

tag_x = marker_w(1,4); 
tag_y = marker_w(2,4); 
r_ = atan2(marker_w(1,1), marker_w(2,1));

%% find x0
w=rodrigues(K(1:3, 1:3));
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K(1,4), K(2,4), K(3,4), r_, tag_x, tag_y];

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'sqp');
x = fmincon(@find_expectedmarker_projection,x0,[],[],[],[],[],[],[], options);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6))

% get refined tb pose
rxy_table = x(7:9) ;
r = rxy_table(1); x=rxy_table(2); y=rxy_table(3);
tb_pose = [ cos(r), -sin(r), 0 , x; ...
            sin(r), cos(r), 0, y ; ...
            0, 0, 1, tb_height; ...
            0 0 0 1]

%% optimization using lsqcurvefit


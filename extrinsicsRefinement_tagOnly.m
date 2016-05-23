addpath(strcat(pwd, '/bkg_code'))

%% load data
load test1_data.mat;
% get projection data from HERB (init estimate) 
load HERB_projectionData; %% A, K, P

load HERB_tableData;

%%

marker_wld_q = tags_file1(1).wld;
marker_w_rot = quatToRotationMatrix(marker_wld_q(4:7))
marker_w_t = marker_wld_q(1:3)';
marker_w = [marker_w_rot marker_w_t; 0 0 0 1];

%expected marker 

table_pose = K * marker_w * tb_kinbody_offset;
r_ = atan2(table_pose(1,1), table_pose(2,1));

%% convert the 3x3 rotation matrix into 3-vector w=[wx wy wz] of the Rodigrues representation
R = K(1:3, 1:3);
theta=acos((trace(R)-1)/2);
w=(theta/(2*sin(theta)))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K(1,4), K(2,4), K(3,4), r_, table_pose(1), table_pose(2) ];

%% optimization 

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
x = fmincon(@find_expectedmarker_projection,x0,[],[]);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

% get refined tb pose
rxy_table = x(7:9) ;
r = rxy_table(1); x=rxy_table(2); y=rxy_table(3);
tb_pose = [ cos(r), -sin(r), 0 , x; ...
            sin(r), cos(r), 0, y ; ...
            0, 0, 1, tb_height; ...
            0 0 0 1];




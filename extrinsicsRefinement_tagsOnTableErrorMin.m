clear
clc
 
%% TRAINING

addpath(strcat(pwd, '/bkg_code'))
addpath(strcat(pwd, '/dataFromTableRotationTest'))

%% LOADING DATA
load HERB_projectionData; %% A, K, P
load('dataFromTableRotationTest/tagTableRotationTest07062016.mat')
load test1_data.mat;
% tags_test
% tags_train

%%  find x

w=rodrigues(K_wc(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4)]

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');
[x, fval, exittag] = fmincon(@minimize_tableTagError,x0,[],[],[],[],[],[],[], options);


% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

K_wc

%% TESTING

load HERB_tableData.mat
zhat_table = [0 0 1]';
z_table = tb_height; 


for i=1:size(tags_test,1) 
    
    
    t = tags_test(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    marker_rot_tabl = (tb_kinbody_offset(1:3, 1:3));
    t_c_ = [inv(marker_rot_tabl)*t_c(1:3, 1:3) t_c(1:3,4) ; 0 0 0 1];
    t_w = K_wc * t_c_;
    
    zhat_tag = t_w(1:3, 3);
    z_tag =  t_w(3, 4);
    
        
    norm_z(i) = norm(z_tag - z_table);
    norm_zhat(i) = norm(dot(zhat_tag,zhat_table));
    
end

norm_z
norm_zhat
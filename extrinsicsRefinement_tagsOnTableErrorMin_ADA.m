clear
clc
 
%% TRAINING

addpath(strcat(pwd, '/bkg_code'))
addpath(strcat(pwd, '/dataFromTableRotationTest'))

%% LOADING DATA
load ADA_Data; %% A, K, P
load('dataFromADA/ADAtags_08062016.mat')


%%  find x

K_wc = T_e2c;

w=rodrigues(K_wc(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4)]

%% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@minimize_tableTagErrorADA,x0,[],[],[],[],[],[],[], options);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1]

K_wc

%% Parameter Refinement

zhat_table = [0 0 1]';
z_table = tag_height; 

% FROM OPTIMIZATION
counter=0;
zhat_tag = [];
for i=1:size(tags_train,1) 
     
    t = tags_train(i, :);
    
    t_c = eye(4,4); 
    t_c(1:3, 4) = t(1:3);
    q = t(4:7);
    t_c(1:3, 1:3) = quatToRotationMatrix(q);
    
    t_w = T_tb2rb*T_rb2e * KK * t_c;
    
    axis_angle = vrrotmat2vec(t_w(1:3, 1:3));
    if(axis_angle(3) <0)
        zhat_tag = -1 .* axis_angle(1:3);
    else
        zhat_tag = axis_angle(1:3);
    end
    
    if(zhat_tag(3) > 0.5)
        counter = counter+1;
    end
    
    z_tag =  t_w(3, 4);
        
    norm_z(i) = norm(z_tag - z_table);
    norm_zhat(i) = dot(zhat_tag,zhat_table);
    
end

%%

z_table = 0.003
% % % % % zhat_table = [0 0 1]'
% % % % % counter = 0
% % % % % for i=1:size(tags_train,1) 
% % % % %     
% % % % %     t = tags_train(i, :);
% % % % %     
% % % % %     t_c = eye(4,4); 
% % % % %     t_c(1:3, 4) = t(1:3);
% % % % %     q = t(4:7);
% % % % %     t_c(1:3, 1:3) = quatToRotationMatrix(q);
% % % % %     
% % % % %     % for ADA
% % % % %     t_w = T_tb2rb*T_rb2e * K_wc * t_c;
% % % % %     qw = rotationMatrixToQ(t_w(1:3, 1:3))';
% % % % %     
% % % % %     axis_angle = vrrotmat2vec(t_w(1:3, 1:3));
% % % % %     if(axis_angle(3) <0)
% % % % %         zhat_tag = -1 .* axis_angle(1:3);
% % % % %     else
% % % % %         zhat_tag = axis_angle(1:3);
% % % % %     end
% % % % %     
% % % % %     if(zhat_tag(3) > 0.5)
% % % % %         z_tag =  t_w(3, 4);
% % % % %         counter = counter + 1;
% % % % %         norm_z(counter) = norm(z_tag - z_table);
% % % % %         norm_zhat(counter) = dot(zhat_tag,zhat_table);
% % % % %     end
% % % % %        
% % % % % end

% disp('****************************')
% disp('Optimized Results')
% norm_z
% norm_zhat
% disp('****************************')
% %% TESTING
% 
% zhat_table = [0 0 1]';
% z_table = tag_height; 
% 
% % FROM OPTIMIZATION
% for i=1:size(tags_test,1) 
%      
%     t = tags_test(i, :);
%     
%     t_c = eye(4,4); 
%     t_c(1:3, 4) = t(1:3);
%     q = t(4:7);
%     t_c(1:3, 1:3) = quatToRotationMatrix(q);
%     
%     t_w = T_tb2rb*T_rb2e * KK * t_c;
%         
%     zhat_tag = t_w(1:3, 3);
%     z_tag =  t_w(3, 4);
%         
%     norm_z(i) = norm(z_tag - z_table);
%     norm_zhat(i) = dot(zhat_tag,zhat_table);
%     
% end
% 
% disp('****************************')
% disp('Optimized Results')
% norm_z
% norm_zhat
% disp('****************************')
% 
% % FROM INITIAL GUESS
% for i=1:size(tags_test,1) 
%      
%     t = tags_test(i, :);
%     
%     t_c = eye(4,4); 
%     t_c(1:3, 4) = t(1:3);
%     q = t(4:7);
%     t_c(1:3, 1:3) = quatToRotationMatrix(q);
%     
%     t_w = T_tb2rb*T_rb2e * K_wc * t_c;
%         
%     zhat_tag = t_w(1:3, 3);
%     z_tag =  t_w(3, 4);
%         
%     norm_z_0(i) = norm(z_tag - z_table);
%     norm_zhat_0(i) = dot(zhat_tag,zhat_table);
%     
% end
% 
% disp('****************************')
% disp('Initial Guess')
% norm_z_0
% norm_zhat_0
% disp('****************************')
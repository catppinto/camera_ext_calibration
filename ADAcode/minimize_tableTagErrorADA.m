function [F] = minimize_tableTagErrorADA(x, tags_train, K_ext0)

%% get projection data from ADA (init estimate) 
load ADA_Data; 
path = '~/cat_workspace/src/ExtrinsicsCalibration/matlab_code/';
% path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/';
% load([path, 'dataFromADA/ADAtags_13062016_trueSizeTag.mat'], 'T_rb2e')
load([path, 'dataFromADA/ADA_depthtests_gt.mat'], 'T_rb2e')

fullcampose = false;
if(length(size(tags_train)) > 2)
    fullcampose = true;
end
    
zhat_table = [0 0 1]';
z_table = tag_height;   

alpha = 10;
gamma = [0.01 0.001]';
psi = 10000;
psi_neg = 1000;
       
A = [567.7969970703125, 0.0, 319.5, 0.0; ...
     0.0, 567.7969970703125, 239.5, 0.0;...
     0.0, 0.0, 1.0, 0.0;...
     0.0, 0.0, 0.0, 1.0];
 
 
 if(fullcampose)
    tag_rotation_offset = eye(4);
else
    tag_rotation_offset = [rotationAroundY(pi) [0;0;0]; 0 0 0 1];
 end


%% Conversion to R,t extrinsiccs 

R_wc = rodrigues(x(1:3));
t_wc=[x(4);x(5);x(6)];

K_wc_est = [R_wc t_wc ; 0 0 0 1];

if(fullcampose)
    s = size(tags_train,3); 
else
    s = size(tags_train,1);
end

%% Load Tag Poses
counter = 0;
for i=1:s
       
    if(fullcampose)
        t = tags_train(:,  :, 1);
        t_c = t;
    else
        t = tags_train(i, :);
        t_c = eye(4,4); 
        t_c(1:3, 4) = t(1:3);
        q = t(4:7);
        t_c(1:3, 1:3) = quatToRotationMatrix(q);
    end
        
    % for ADA
    t_w = T_tb2rb * T_rb2e * K_wc_est * t_c;
    t_w = tag_rotation_offset * t_w;
    qw = rotationMatrixToQ(t_w(1:3, 1:3))';
    t_w_estimates(:, :, i) = t_w;
    
    zhat_tag = t_w(1:3, 3);
    
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_z(counter) = norm(z_tag - z_table);
    
    if(z_tag < 0)
        negativeZ(counter) = abs(z_tag);
    else
        negativeZ(counter) = 0; 
    end;
    
    CosTheta = dot(zhat_tag,zhat_table)/(norm(zhat_tag)*norm(zhat_table));
    theta_error(counter) = acos(CosTheta);
       
end
    
r_ext0 = K_ext0(1:3); 
norm_rExt = norm(r_ext0 - x(1:3));
t_ext0 = K_ext0(4:6); 
norm_tExt = norm(t_ext0 - x(4:6));
norm_ext = [norm_rExt ; norm_tExt];


%% cost function 
F = sum(psi*norm_z)+ psi_neg*sum(negativeZ) + sum(alpha*theta_error); % + sum(gamma.*norm_ext); % + delta*reprojection_error; 
% delta = 1;
% F = sum(psi*norm_z)+ psi_neg*sum(negativeZ) + sum(alpha*theta_error) + delta*reprojection_error; 



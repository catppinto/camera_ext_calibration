function [F] = minimize_tableTagErrorADA(x, tags_train, K_ext0, transforms)

%% get projection data from ADA (init estimate) 

T_w2rb = transforms.T_w2rb;
T_rb2e = transforms.T_rb2e;
tag_height = transforms.tag_height;

fullcampose = false;
if(length(size(tags_train)) > 2)
    fullcampose = true;
end
    
zhat_table = [0 0 -1]';
z_table = tag_height;   

alpha = 0.01;
gamma = [1 1 0.01 0.1 0.1 0.001]';
psi = 10;
 
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
    t_w = T_w2rb * T_rb2e * K_wc_est * t_c;
    
    zhat_tag = t_w(1:3, 3);
    
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_z(counter) = norm(z_tag - z_table);
       
    CosTheta = dot(zhat_tag,zhat_table)/(norm(zhat_tag)*norm(zhat_table));
    theta_error(counter) = abs(acos(CosTheta));
       
end

x(1)
K_ext0(1)  
norm_rExt = [ ( x(1) - K_ext0(1)).^2 ; 
              ( x(2) - K_ext0(2)).^2 ;
             abs( (K_ext0(3) - x(3)))];
t_ext0 = K_ext0(4:6); 
norm_tExt = [abs( (t_ext0(1) - x(4))) ; 
             abs( (t_ext0(2) - x(5))) ;
             abs( (t_ext0(3) - x(6)))];
norm_ext = [norm_rExt ; norm_tExt];

%% cost function 

norm_z
theta_error
norm_ext
 
F = sum(psi*norm_z) + sum(alpha*theta_error) + sum(gamma.*norm_ext) ;



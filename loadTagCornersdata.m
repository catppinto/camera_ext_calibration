function [tagcorners, halftagsize] = loadTagCornersdata(filename) 

load (filename) ;
load HERB_projectionData; %% A, K, P

% load tagcenter in uv coordinates (image frame)
uv_c1 = tags_file1(1).cam_corners(1,1:3);
uv_c2 = tags_file1(1).cam_corners(2,1:3);
uv_c3 = tags_file1(1).cam_corners(3,1:3);
uv_c4 = tags_file1(1).cam_corners(4,1:3);

% load tagcenter in cam coordinates
tagcenter_Xc = eye(4,4);
tag_pose_camFrame = tags_file1(1).wld_cam;
tagcenter_Xc(1:3, 1:3) = quatToRotationMatrix(tag_pose_camFrame(4:7));
tagcenter_Xc(1:3, 4) = tag_pose_camFrame(1:3);

% convert tagcenter_Xc to tagcenter_Xw
tagcenter_Xw = K_wc * tagcenter_Xc;

%find wld coordinates of corners
halftagsize = tabletag_size/2;
w_m  = [tagcenter_Xw(1:3,4)];
w_rot = tagcenter_Xw(1:3, 1:3);

tagcorners.uv.c1 = uv_c1;
tagcorners.uv.c2 = uv_c2;
tagcorners.uv.c3 = uv_c3;
tagcorners.uv.c4 = uv_c4;
tagcorners.xyz_w = w_m;
tagcorners.rot_w = w_rot;


% marker_rot_tabl = (tb_kinbody_offset(1:3, 1:3));
% %Xc_ = [inv(marker_rot_tabl)*marker_wc_rot marker_wc_t ; 0 0 0 1]
% %Xw = K_wc * Xc_ 
% Xw = K_wc * Xc
clear 
clc


path = '/home/cat/Documents/tese/FEUP_camera_calibration/DATA/';
addpath(path)
fid = fopen('/home/cat/catkin_ws/src/biovision/cat_robotcam_calibration/markerarray_log.txt','r')


C = textscan(fid, '%s','Delimiter','');
fclose(fid)
C = C{:};

[x,y]=size(C);

index = strfind(C,'position:');
index = find(~cellfun(@isempty,index));

ntags = length(index);
tag = [];
for i = 1:ntags
    
    pos_x = C{index(i)+1};
    pos_x = str2double(pos_x(regexp(pos_x,[': '])+1:end));
    pos_y = C{index(i)+2};
    pos_y = str2double(pos_y(regexp(pos_y,[': '])+1:end));
    pos_z = C{index(i)+3};
    pos_z = str2double(pos_z(regexp(pos_z,[': '])+1:end));
    ori_x = C{index(i)+5};
    ori_x = str2double(ori_x(regexp(ori_x,[': '])+1:end));
    ori_y = C{index(i)+6};
    ori_y = str2double(ori_y(regexp(ori_y,[': '])+1:end));
    ori_z = C{index(i)+7};
    ori_z = str2double(ori_z(regexp(ori_z,[': '])+1:end));
    ori_w = C{index(i)+8};
    ori_w = str2double(ori_w(regexp(ori_w,[': '])+1:end));
        
    tag(i, :) = [pos_x pos_y pos_z ori_x ori_y ori_z ori_w];
    
end  

index1 = strfind(C,'corners2d: ');
index1 = find(~cellfun(@isempty,index1));

for i = 1:ntags
    
    corner1_x = C{index1(i)+2};
    corner1_x = str2double(corner1_x(regexp(corner1_x,[': '])+1:end));
    corner1_y = C{index1(i)+3};
    corner1_y = str2double(corner1_y(regexp(corner1_y,[': '])+1:end));
    
    corner2_x = C{index1(i)+6};
    corner2_x = str2double(corner2_x(regexp(corner2_x,[': '])+1:end));
    corner2_y = C{index1(i)+7};
    corner2_y = str2double(corner2_y(regexp(corner2_y,[': '])+1:end));
    
    corner3_x = C{index1(i)+10};
    corner3_x = str2double(corner3_x(regexp(corner3_x,[': '])+1:end));
    corner3_y = C{index1(i)+11};
    corner3_y = str2double(corner3_y(regexp(corner3_y,[': '])+1:end));
    
    corner4_x = C{index1(i)+14};
    corner4_x = str2double(corner4_x(regexp(corner4_x,[': '])+1:end));
    corner4_y = C{index1(i)+15};
    corner4_y = str2double(corner4_y(regexp(corner4_y,[': '])+1:end));
        
    corners_1(i, :) = [corner1_x corner1_y 1.0];
    corners_2(i, :) = [corner2_x corner2_y 1.0];
    corners_3(i, :) = [corner3_x corner3_y 1.0];
    corners_4(i, :) = [corner4_x corner4_y 1.0];  
end  

ids = [0 1 2 11 4 5 6 7];
wld_points = [0.7005 0.024 0; ...
	0.4 0.024 0; ...
	0.2495 0.024 0; ...
	0.402 0.19 0; ...
	0.449 -0.004 0.16627; ...
	0.449 -0.004 0.11877; ...
	0.449 -0.004 0.07127; ...
	0.449 -0.004 0.02377];

numtags = length(index);
calibration_tags.cam_points = tag; 
calibration_tags.id = ids;
calibration_tags.wld_points = wld_points;
calibration_tags.corners.c1 = corners_1;
calibration_tags.corners.c2 = corners_2;
calibration_tags.corners.c3 = corners_3;
calibration_tags.corners.c4 = corners_4;

save([path,'/FEUPtags.mat'], 'calibration_tags')

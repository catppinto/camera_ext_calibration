clear 
clc

%addpath('/home/cat/CMU_TESTS_DATA/ADAtags_data')
%fid = fopen('adatags_08062016.txt','r')

path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/dataFromADA/';
addpath(path)
fid = fopen('adatags_13062016.txt','r')


C = textscan(fid, '%s','Delimiter','');
fclose(fid)
C = C{:};

[x,y]=size(C);

index = strfind(C,'position:');
index = find(~cellfun(@isempty,index));


tag = [];
for i = 1:length(index)
    
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

index1 = strfind(C,'corners2d[]');
index1 = find(~cellfun(@isempty,index1));

for i = 1:length(index1)
    
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

numtags = length(index);
tags_train.pose = tag(1:round(numtags/2), :); 
tags_test.pose = tag(round(numtags/2+1):end, :);
tags_train.corners.c1 = corners_1(1:round(numtags/2), :);
tags_train.corners.c2 = corners_2(1:round(numtags/2), :);
tags_train.corners.c3 = corners_3(1:round(numtags/2), :);
tags_train.corners.c4 = corners_4(1:round(numtags/2), :);
tags_test.corners.c1 = corners_1(round(numtags/2+1):end, :);
tags_test.corners.c2 = corners_2(round(numtags/2+1):end, :);
tags_test.corners.c3 = corners_3(round(numtags/2+1):end, :);
tags_test.corners.c4 = corners_4(round(numtags/2+1):end, :);

save([path,'/ADAtags_13062016.mat'], 'tags_train', 'tags_test')

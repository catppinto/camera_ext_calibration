%addpath('/home/cat/CMU_TESTS_DATA/ADAtags_data')
%fid = fopen('adatags_08062016.txt','r')

path = '/home/apirespi/cat_workspace/src/tabletop_perception_tools/pcd_test_files/ADA_sameRotationTags_11062016';
addpath(path)
fid = fopen('adatags_11062016.txt','r')


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

tags_train = tag(1:round(i/2), :); 
tags_test = tag(round(i/2+1):end, :);

save([path,'/ADAtags_11062016.mat'], 'tags_train', 'tags_test')

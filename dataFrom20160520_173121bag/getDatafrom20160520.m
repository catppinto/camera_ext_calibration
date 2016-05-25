%% tag calibration data from tag_detection file 
% /home/cat/Documents/CMU_Herb/camera_ext_calibration/dataFrom20160520_173121bag

function [cam_points, wld_points] = getDatafrom20160520()
c1 = [437.5, 345.5, 1.0];
c2 = [410.5, 345.5, 1.0];
c3 = [407.5, 362.5, 1.0];
c4 = [436.5, 362.5, 1.0];
c14 = [(c1(1)+c4(1))/2 (c1(2)+c4(2))/2 1];
c23 = [(c2(1)+c3(1))/2 (c2(2)+c3(2))/2 1];
m = [423 354 1];

tag_size= 0.0480000004172;
ts = 0.0480000004172/2;

w_m =  [-0.0753773175125, 0.137904712863, 0.916102764935 1]

w_c1 = [w_m(1)+ts w_m(2)+ts w_m(3) 1];
w_c2 = [w_m(1)-ts w_m(2)+ts w_m(3) 1];
w_c3 = [w_m(1)-ts w_m(2)-ts w_m(3) 1];
w_c4 = [w_m(1)+ts w_m(2)-ts w_m(3) 1];
w_c23 = [w_m(1)-ts w_m(2) w_m(3) 1];
w_c14 = [w_m(1)+ts w_m(2) w_m(3) 1];

% figure;
% scatter3(w_m(1), w_m(2), w_m(3), 'b', 'filled')
% hold on 
% scatter3(w_c1(1), w_c1(2), w_c1(3), 'b')
% hold on 
% scatter3(w_c2(1), w_c2(2), w_c2(3), 'g')
% hold on
% scatter3(w_c3(1), w_c3(2), w_c3(3), 'y')
% hold on
% scatter3(w_c4(1), w_c4(2), w_c4(3), 'm')
% hold on
% scatter3(w_c14(1), w_c14(2), w_c14(3), 'p')
% hold on
% scatter3(w_c23(1), w_c23(2), w_c23(3), 'p')
% hold on

cam_points = [c1; c2; c3; c4; c14; c23; m];

wld_points = [w_c1; w_c2; w_c3; w_c4; w_c14; w_c23; w_m];
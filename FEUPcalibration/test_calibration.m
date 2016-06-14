% close all
% clc
% clear all

addpath(genpath('/home/cat/Documents/tese/HOPE_matlab'));
n= 4;

w = [0.7005 0.024 0; ... 
   0.4 0.024 0 ; ...
   0.2495 0.024 0; ...
   0.402 0.19 0; ...
   0.449 -0.004 0.16627; ...
   0.449 -0.004 0.11877; ...
   0.449 -0.004 0.07127; ...
   0.449 -0.004 0.02377];

ctag = [-0.322713 0.0119508 0.671776; ...
    0.0577349 0.0274597 0.670989; ...
    0.243119 0.0353418 0.669769; ...
    0.053183 0.158643 0.511146; ...
    0.0110721 -0.161089 0.57441; ...
    0.00809416 -0.114525 0.611413; ...
    0.00515552 -0.0678623 0.646763; ...
    0.00249965 -0.0216926 0.685821];
ctag_h = [ctag ones(8,1)]';

% scatter3(w(:, 1), w(:, 2), w(:, 3), 'r');
% hold on;
% scatter3(ctag(:, 1), ctag(:, 2), ctag(:, 3), 'g');

T = performCalibration(w, ctag); 

%% point in world frame
[w_p] = getWorldPointFromCameraPoint(ctag, T);
w_test = w_p(n, :); 

fprintf('w_test : %3.7f, %3.7f, %3.7f \n', ...
        w_test(1), w_test(2), w_test(3));

% figure;
% hold on
% scatter3(w(:,1), w(:,2), w(:,3), 'r');
% hold on
% scatter3(w_p(:,1), w_p(:, 2), w_p(:, 3), 'b');
% a = ww * 1.25;
% hold on
% scatter3(a(1,:), a(2,:), a(3,:), 'p');
   
%% point in manipulator frame
T_worldToManip = [ eye(3), [0.0135 0.0135 0.023]'; 0 0 0 1];
w_test_manip = inv(T_worldToManip) * w_test'; 
w_test_manip = w_test_manip(1:3, 1);

fprintf('w_test_manip : %3.7f, %3.7f, %3.7f \n', ...
        w_test_manip(1), w_test_manip(2), w_test_manip(3));
    

%% perform IK
pos = w_test_manip;
rot = rotationAroundX(pi);
  
[theta_vector, s] = compute_ik(pos, rot);
if(s)
    
    theta_vector(3) = theta_vector(3) + pi/4; 

    [is_sing_J11, is_sing_J22, J, J11, J22] = find_singularities(theta_vector);

    if ( is_sing_J11 )
        disp('Has a singularity in the first part of the manipulator'); 
    elseif (is_sing_J22)
        disp('Has a singularity in the wrist joints'); 
    else

        [eef_xyz, eef_rot] = forward_kin_moreira(theta_vector);  

        fprintf('GT xyz : %3.7f, %3.7f, %3.7f \n', ...
            pos(1), pos(2), pos(3));

        fprintf('PD xyz : %3.7f, %3.7f, %3.7f \n', ...
            eef_xyz(1), eef_xyz(2), eef_xyz(3));

        fprintf('PD joints : %3.7f, %3.7f, %3.7f, %3.7f, %3.7f, %3.7f \n', ...
            theta_vector(1), theta_vector(2), theta_vector(3), ...
            theta_vector(4), theta_vector(5), theta_vector(6));
    end
else
    fprintf('Not reachable point \n');
end


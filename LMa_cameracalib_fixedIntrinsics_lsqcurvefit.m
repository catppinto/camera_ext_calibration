% uv = c.cam_points;
% w = c.wld_points;
% 
% uv_x = uv(:, 1)./uv(:, 3);
% uv_y = uv(:, 2)./uv(:, 3);
% 
% % project wld points
% w_p = (c.proj_matrix_computed * w')';
% 
% w_p_x = w_p(:, 1)./w_p(:, 3);
% w_p_y = w_p(:, 2)./w_p(:, 3);
% 
% d_x = uv_x - w_p_x
% d_y = uv_y - w_p_y


clc
clear 

c = CamExtrinsicsCalibrationClass();

n = 100;
wx = 0 + (2-0).*rand(n,1); 
wy = 0 + (2-0).*rand(n,1); 
wz = 0 + (5-0).*rand(n,1);

wld_points = [wx wy wz ones(n, 1)];

c = c.CalculateExtrinsics(wld_points, true);

%% set initial guesses
% t0 = c.t_ext;
% R0 = c.R_ext;

t0 = c.extrinsic_camera_matrix(1:3, 4);
R0 = c.extrinsic_camera_matrix(1:3, 1:3);



%% Using matlab function lsqcurvefit

xdata = c.wld_points';
uv = c.cam_points';
ydata(1, :) = uv(1, :)./uv(3, :);
ydata(2, :) = uv(2, :)./uv(3, :);

% R0, t0 = initial extrinsic parameters obtained from the solution in Section 3
tx=t0(1);
ty=t0(2);
tz=t0(3);

% convert the 3x3 rotation matrix into 3-vector w=[wx wy wz] of the Rodigrues representation
R=R0;
theta=acos((trace(R)-1)/2);
w=(theta/(2*sin(theta)))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, tx, ty, tz];

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Display', 'iter');
lb = [];
ub = [];
[x1,resnorm,residual,exitflag,output]  = lsqcurvefit(@myfun,x0,xdata,ydata, lb,ub,options);
[R_est1, t_est1] = xToRt(x1);


options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Jacobian','on', 'Display', 'iter');
lb = [];
ub = [];
[x2,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(@myfun,x0,xdata,ydata, lb,ub,options);
[R_est2, t_est2] = xToRt(x2);


%%  displays
K = c.intrinsic_camera_matrix;
disp(' Ground Truth ') 
c.extrinsic_camera_matrix

disp(' Initial Guess ') 
m = [c.R_ext c.t_ext; 0 0 0 1]

disp(' LM no Jacobian ') 
m = [R_est1 t_est1; 0 0 0 1]
p = K * m;
error_1 = computeReprojectionError( c.wld_points, c.proj_matrix, p);

disp(' LM Jacobian on ') 
m = [R_est2 t_est2; 0 0 0 1]
p = K * m;
error_2 = computeReprojectionError( c.wld_points, c.proj_matrix, p);

fprintf('Error no JACOBIAN : %f \n Error w JACOBIAN : %f \n', error_1, error_2);  

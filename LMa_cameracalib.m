%%% https://engineering.purdue.edu/kak/computervision/ECE661/HW5_LM_handout.pdf



%% LMa for Camera Calibration Parameters Refinement

% % au = alpha, bu = beta, sk = gamma
% syms u0 v0 au av real
% syms tx ty tz wx wy wz real
% syms X Y Z um vm real
% 
% % intrinsics
% K = [ au 0 u0 0; ...
%     0 av v0 0; ...
%     0 0 1 0;
%     0 0 0 1];
% 
% 
% % expression for rotation matrix from rodrigues formula
% theta = sqrt(wx^2 + wy^2 + wz^2);
% omega = [   0 -wz wy; ...
%             wz 0 -wx; ...
%             -wy wx 0];
% R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
% 
% % translation vector
% t = [tx; ty; tz];
% 
% % perspective projection of the model point (X,Y)
% p_matrix = K*[R t; 0 0 0 1];
% uvs=p_matrix*[X; Y; Z; 1];
% u=uvs(1)/uvs(3);
% v=uvs(2)/uvs(3);
% 
% % calculate the geometric distance in x and y direction
% % um,vm =the x and y positions of an extracted corner
% % u,v = the x and y positions of the projection of the corresponding model point
% dx=um-u;
% dy=vm-v;
% 
% % Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
% Jx=jacobian(dx,[au,av,u0,v0,wx, wy, wz, tx, ty, tz]);
% Jy=jacobian(dy,[au,av,u0,v0,wx, wy, wz, tx, ty, tz]);

%% 
clc
clear 

c = CamExtrinsicsCalibrationClass();

n = 100;
wx = 0 + (2-0).*rand(n,1); 
wy = 0 + (2-0).*rand(n,1); 
wz = 0 + (2-0).*rand(n,1);

wld_points = [wx wy wz ones(n, 1)];

c = c.CalculateExtrinsics(wld_points, true);

% set initial guesses
t0 = c.t_ext;
R0 = c.R_ext;
K0 = c.intrinsic_camera_matrix;
au = c.intrinsic_camera_matrix(1,1);
u0 = c.intrinsic_camera_matrix(1,3);
av = c.intrinsic_camera_matrix(2,2);
v0 = c.intrinsic_camera_matrix(2,3);
xe = (c.cam_points(:, 1:3))';
Xw = (c.wld_points(:, 1:3))';

%% LM aplication

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


iters=100; 				% set the number of iterations for the LM algorithm
lamda=0.1; 				% initial the value of damping factor
updateJ=1;
Ndata=2*size(c.cam_points, 1);
Nparams=10;
for it=1:iters
    if updateJ==1

        % create the intrinsic parameter matrix
        K=[ au 0  u0 0;
            0  av v0 0;
            0  0   1 0;
            0  0   0 1];

        % convert the 3-vector [wx wy wz] of the Rodigrues representation
        % into the 3x3 rotation matrix
        theta2=wx^2+wy^2+wz^2;
        theta=sqrt(theta2);
        omega=[0 -wz wy;
            wz 0  -wx;
            -wy wx 0;];
        R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
        t=[tx;ty;tz];

        % Evaluate the Jacobian at the current parameter values
        % and the values of geometric distance
        J=zeros(Ndata,Nparams);
        d=zeros(Ndata,1);
        for i=1:size(xe,2)
            X=Xw(1,i); 		% (X,Y) are the coordinates of the ith model point
            Y=Xw(2,i);
            Z=Xw(3,i);
            J(2*(i-1)+1,:)= [ -(tx - Y*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + X*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)), 0, -1, 0, (X*(au*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - u0*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - Y*(au*((wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + u0*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(u0*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - au*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) + ((X*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(au*tx - Y*(au*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - au*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(au*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, (X*(au*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - u0*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(u0*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - au*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) - Y*(au*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + u0*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) + ((Y*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + X*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(au*tx - Y*(au*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - au*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(au*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, (Z*(u0*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - au*((wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - Y*(u0*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + au*((wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wz^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + X*(au*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - u0*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) + ((X*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2))*(au*tx - Y*(au*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - au*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(au*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, -au/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)), 0, (au*tx - Y*(au*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - u0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*u0 - X*(u0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - au*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wy^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) + Z*(au*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + u0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2 - u0/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))];
            J(2*(i-1)+2,:)= [ 0, -(ty + X*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - Z*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)), 0, -1, (Y*(av*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - v0*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(v0*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - av*((wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) - X*(av*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx^2*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + v0*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) + ((X*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(av*ty + X*(av*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + av*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(av*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, (Y*(av*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - v0*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - X*(av*((wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wx*wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + v0*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Z*(v0*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - av*((wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) + ((Y*((wx*wy*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy^2*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy^2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + X*((wy^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2)))*(av*ty + X*(av*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + av*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(av*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, (Z*(v0*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - av*((wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)) - X*(v0*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) + av*(sin((wx^2 + wy^2 + wz^2)^(1/2))/(wx^2 + wy^2 + wz^2)^(1/2) + (wz^2*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2 + (wx*wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2))) + Y*(av*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^2 - (2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2)^(3/2)) - v0*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)) + ((X*((wy*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) - (wy*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wx*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) + (wx*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wx*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2) - Z*((2*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^2 + (wz*sin((wx^2 + wy^2 + wz^2)^(1/2))*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2)^(3/2)) + Y*((wx*wz*cos((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2) - (wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2) - (wx*wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (wy*wz^2*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(3/2) + (2*wy*wz^2*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)^2))*(av*ty + X*(av*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + av*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(av*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2, 0, -av/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)), (av*ty + X*(av*((wz*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wx*wy*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2))) + tz*v0 + Y*(v0*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + av*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wz^2))/(wx^2 + wy^2 + wz^2) + 1)) - Z*(av*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) - v0*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1)))/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))^2 - v0/(tz - X*((wy*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) + (wx*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Y*((wx*sin((wx^2 + wy^2 + wz^2)^(1/2)))/(wx^2 + wy^2 + wz^2)^(1/2) - (wy*wz*(cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1))/(wx^2 + wy^2 + wz^2)) + Z*(((cos((wx^2 + wy^2 + wz^2)^(1/2)) - 1)*(wx^2 + wy^2))/(wx^2 + wy^2 + wz^2) + 1))];

            % perspective projection of the model point
            p_matrix = K*[R t; 0 0 0 1];
            uvs=p_matrix*[X; Y; Z; 1];
            up=uvs(1)/uvs(3);
            vp=uvs(2)/uvs(3);
            % compute the geometric distance in x & y directions
            % xe(1,i), xe(2,i) = the the x and y positions of the ith extracted corner.
            % up,vp = the x and y positions of the projection of the corresponding model point
            d(2*(i-1)+1,1)=xe(1,i)-up;
            d(2*(i-1)+2,1)=xe(2,i)-vp;
        end

        % compute the approximated Hessian matrix
        H=J'*J;
        if it==1 % the first iteration : compute the initial total error
            e=dot(d,d);
            disp(e);
        end
    end

    % Apply the damping factor to the Hessian matrix
    H_lm=H+(lamda*eye(Nparams,Nparams));

    % Compute the updated parameters
    dp=-inv(H_lm)*(J'*d(:));
    au_lm=au+dp(1);
    av_lm=av+dp(2);
    u0_lm=u0+dp(3);
    v0_lm=v0+dp(4);
    wx_lm=wx+dp(5);
    wy_lm=wy+dp(6);
    wz_lm=wz+dp(7);
    tx_lm=tx+dp(8);
    ty_lm=ty+dp(9);
    tz_lm=tz+dp(10);

    % Evaluate the total geometric distance at the updated parameters
    K=[au_lm 0     u0_lm 0;
        0    av_lm v0_lm 0;
        0    0       1   0;
        0    0       0   1];
    omega=[0 -wz_lm wy_lm;
        wz_lm 0  -wx_lm;
        -wy_lm wx_lm 0;];
    theta2=wx_lm^2+wy_lm^2+wz_lm^2;
    theta=sqrt(theta2);
    R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
    t=[tx_lm;ty_lm;tz_lm];
    d_lm=zeros(Ndata,1);
    for i=1:size(xe,2)
        X=Xw(1,i);
        Y=Xw(2,i);
        Z=Xw(3,i);
        p_matrix = K*[R t; 0 0 0 1];
        uvs=p_matrix*[X; Y; Z; 1];
        up=uvs(1)/uvs(3);
        vp=uvs(2)/uvs(3);
        d_lm(2*(i-1)+1,1)=xe(1,i)-up;
        d_lm(2*(i-1)+1,1)=xe(2,i)-vp;
    end

    % e_lm is the error between the image coordinates and projective coordinates using
    % the updated parameters
    e_lm=dot(d_lm,d_lm);

    %If the total geometric distance of the updated parameters is less than the previous one
    % then makes the updated parameters to be the current parameters
    % and decreases the value of the damping factor
    if e_lm<e
        lamda=lamda/10;
        au=au_lm;
        av=av_lm;
        u0=u0_lm;
        v0=v0_lm;
        wx=wx_lm;
        wy=wy_lm;
        wz=wz_lm;
        tx=tx_lm;
        ty=ty_lm;
        tz=tz_lm;
        e=e_lm;
        disp(e);
        updateJ=1;
    else
        % otherwise increase the value of the damping factor and try again
        updateJ=0;
        lamda=lamda*10;
    end
end


%%  displays

disp(' Ground Truth ') 
c.extrinsic_camera_matrix

disp(' Initial Guess ') 
c.R_ext
c.t_ext

disp(' LM Results ') 
R
t


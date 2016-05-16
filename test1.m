     clear
     clc
     
cam_points = [ ...
    0.1    0.2    1.0000    1.0000 ; ...
    0.3    0.4    1.0000    1.0000 ; ...
    0.7    0.5    1.0000    1.0000 ; ...
    0.5    0.4    1.0000    1.0000 ; ...
    0.2    0.3    1.0000    1.0000 ; ...
    0.3    0.2    1.0000    1.0000 ; ...
    0.4    0.1    1.0000    1.0000 ; ...
    0.5    0.2    1.0000    1.0000 ; ...
    0.1    0.2    1.0000    1.0000 ; ...
    0.2    0.3    1.0000    1.0000];

% n = 100;
% cam_u = 0 + (1-0).*rand(n,1); 
% cam_v = 0 + (1-0).*rand(n,1); 
% cam_z = ones(n, 1); 
% 
% cam_points = [cam_u cam_v cam_z ones(n, 1)];
% 
% p = [ ...
%      1    0    0   1; ...
%      0    1    0   1; ...
%      0    0    1   1; ...
%      0    0    0   1; ...
%      ];
% 
%  wld_points = (pinv(p)*cam_points')';
% wld_points = round(wld_points, 5);


n = 10;
wx = 0 + (2-0).*rand(n,1); 
wy = 0 + (2-0).*rand(n,1); 
wz = 0 + (2-0).*rand(n,1);

wld_points = [wx wy wz ones(n, 1)];


% p_r = rotationAroundZ(1.2)*rotationAroundY(1.3)*rotationAroundX(-0.5);
% p_t= [5 -3 2]';
% 
% p = [ ...
%      p_r p_t;
%      0    0    0   1; ...
%      ];

ext_t =  [-0.0581, -0.0993, 1.267]';
ext_r = [1.792e-04, -0.357, 0.934;  ...
         -0.1, -6.0504e-04, -3.94e-05; ...
         5.79244942e-04, -9.34119069e-01, -3.56961103e-01];
extrinsic_camera_matrix = [ext_r ext_t; 0 0 0 1];

intrinsic_camera_matrix = ...
    [530, 0.0, 470, 0.0 ; ...
    0.0, 530, 270 , 0.0 ; ...
    0.0, 0.0, 1.0, 0.0;
    0.0, 0.0, 0.0, 1.0];

p = intrinsic_camera_matrix * extrinsic_camera_matrix;


%find camera points
cam_points = (p*wld_points')';

% G matrix
G_matrix = [];
for i=1:size(cam_points, 1)
    uv_cam = cam_points(i,:);
    xyz_world = wld_points(i,:);

    u_cam = uv_cam(1) / uv_cam(3);
    v_cam = uv_cam(2) / uv_cam(3);

    Gentry = [xyz_world(1), xyz_world(2), xyz_world(3), 1, ...
        0, 0, 0, 0, ...
        -u_cam*xyz_world(1), -u_cam*xyz_world(2), -u_cam*xyz_world(3), -u_cam; ...
        0, 0, 0, 0, ...
        xyz_world(1), xyz_world(2), xyz_world(3), 1, ...
        -v_cam*xyz_world(1), -v_cam*xyz_world(2), -v_cam*xyz_world(3), -v_cam];

    G_matrix = [G_matrix; Gentry];
end
            

% compute P matrix : the eigenvector of GtG with smallest
% eigenvalue
[U, S, V] = svd(G_matrix);
eigenv_gtg = V(:, end);

pp = [ eigenv_gtg(1) eigenv_gtg(2) eigenv_gtg(3) eigenv_gtg(4) ; ...
                eigenv_gtg(5) eigenv_gtg(6) eigenv_gtg(7) eigenv_gtg(8) ; ...
                eigenv_gtg(9) eigenv_gtg(10) eigenv_gtg(11) eigenv_gtg(12) ; ...
                0                  0                  0                  1                  ];

pp./(sqrt(pp(3,1)^2+pp(3,2)^2+pp(3,3)^2))

p
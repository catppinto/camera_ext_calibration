clc
clear 


c = CamExtrinsicsCalibrationClass();

% cam_points = [ ...
%     %tag
%     530 192 1; ...
%     515 193 1; ...
%     529 191 1; ...
%     523 193 1; ...
%     531 196 1; ...
%     528 192 1; ...
%     % glass
%     457 138 1; ...
%     466 139 1; ...
%     426 147 1; ...
%     456 137 1; ...
%     465 138 1; ...
%     427 146 1; ...
%     ];

% temp =cam_points(:, 1);
% cam_points(:, 1) = cam_points(:, 2);
% cam_points(:, 1) = temp;

cam_u = 0 + (1-0).*rand(100,1); 
cam_v = 0 + (1-0).*rand(100,1); 
cam_z = ones(100, 1); 
cam_points = [cam_u cam_v cam_z];

c = c.CalculateExtrinsics(cam_points);

%%

            obj = c

            m = obj.proj_matrix_computed;
            r = zeros(3,3);
            t = zeros(3,1); 
            a = eye(3,3);

            gamma = sqrt(m(3,1)^2 + m(3,2)^2 + m(3,3)^2)

            B = m(1:3, 1:3); 
            B = B/gamma;
            b = m(1:3, 4); 

            k = B*B';
            k2 = c.intrinsic_camera_matrix* c.intrinsic_camera_matrix';
            
            %u0
            a(1,3) = k(1,3); 
            %v0
            a(2,3) = k(2,3);
            %beta 
            a(2,2) = sqrt(k(2,2)-a(2,3)^2);
            %alpha
            a(1,1) = sqrt(k(1,1)-a(1,3)^2);

            r = inv(a) * B 
            [U, ~, V] = svd(r);
            d = sign(det(V*U'));
            r = [1 0 0; 0 1 0; 0 0 d]*r;

            t = inv(a) * b
            
            intr = c.intrinsic_camera_matrix(1:3, 1:3);
            t2 = inv(intr) * b
                     
            
            obj.R_noint = r;
            obj.t_noint = t;
            
            
            c = obj;
            









clc
clear 


c = CamExtrinsicsCalibrationClass();

% cam_points = [559 177 1; ...
%               537 171 1; ...
%               574 178 1; ...
%               422 180 1; ...
%               488 187 1; ...
%               564 180 1; ...
%               551 174 1; ...
%               560 186 1; ...
%               553 173 1; ...
%               595 181 1; ...
%               576 179 1; ...
%               524 195 1; ...
%               522 192 1; ...
%               531 192 1; ...
%               534 194 1; ...
%               605 186 1];


cam_points = [ ...
    %tag
    530 192 1; ...
    515 193 1; ...
    529 191 1; ...
    % glass
    457 138 1; ...
    466 139 1; ...
    426 147 1; ...
    ];

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
            %u0
            a(1,3) = k(1,3); 
            %v0
            a(2,3) = k(2,3);
            %beta 
            a(2,2) = sqrt(k(2,2)-a(2,3)^2);
            %alpha
            a(1,1) = sqrt(k(1,1)-a(1,3)^2);
            a

            r = inv(a) * B 
            [U, ~, V] = svd(r);
            d = sign(det(V*U'));
            r_ = [1 0 0; 0 1 0; 0 0 d]*r;

            t = inv(a) * b
            
            tz = m(3,4); 
%             tz= -tz;
            tx = (m(1,4) - a(1,3)*tz)/a(1,1); 
            ty = (m(2,4) - a(2,3)*tz)/a(2,2);
            ta = [tx; ty; tz]
            
            intr = c.intrinsic_camera_matrix;
            tz = m(3,4); 
%             tz= -tz;
            tx = (m(1,4) - intr(1,3)*tz)/intr(1,1);
            ty = (m(2,4) - intr(2,3)*tz)/intr(2,2);
            t2 = [tx; ty; tz]
            
            
            obj.R_noint = r;
            obj.t_noint = t;
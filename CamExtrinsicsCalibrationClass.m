classdef CamExtrinsicsCalibrationClass
    properties
        % initial camera matrix
        intrinsic_camera_matrix
        % extrinsics camera matrix (initial guess)
        extrinsic_camera_matrix
        % proj matrix initial guess
        proj_matrix
        
        cam_points
        wld_points
        
        % G matrix
        G_matrix
        %
        proj_matrix_computed
       
        R_int
        R_int_n
        t_int
        R_noint
        R_noint_n
        t_noint
        
    end
    methods
        function obj = CamExtrinsicsCalibrationClass()

            ext_t =  [-0.058121267086213846, -0.09935999999999337, 1.2668267743218642]';
            ext_r = [1.79185474e-04, -3.56961060e-01, 9.34119248e-01;  ...
                     -9.99999816e-01, -6.05046094e-04, -3.93873210e-05; ...
                     5.79244942e-04, -9.34119069e-01, -3.56961103e-01];
            obj.extrinsic_camera_matrix = [ext_r ext_t; 0 0 0 1];
            
            obj.intrinsic_camera_matrix = ...
                [529.2945040622658, 0.0, 466.96044871160075, 0.0 ; ...
                0.0, 531.2834529497384, 273.2593671723483 , 0.0 ; ...
                0.0, 0.0, 1.0, 0.0;
                0.0, 0.0, 0.0, 1.0];
            
            obj.proj_matrix = obj.intrinsic_camera_matrix * obj.extrinsic_camera_matrix;

        end
        
        function obj = CalculateExtrinsics(obj, cam_points)
            
            obj.cam_points = cam_points;
            
            obj.wld_points = obj.uvToXYZwolrd(cam_points);
            
            obj = obj.constructPmatrix();
            
        end
        
        function xyz_points = uvToXYZwolrd(obj, uv_points)
            
            uv_points = [uv_points ones(size(uv_points, 1),1)];
            
            xyz_points = (inv(obj.proj_matrix) *  uv_points')'; 
            xyz_points = xyz_points(:, 1:3);
            
            for i=1:1:size(uv_points,1)
               fprintf('uv  : %3.4f, %3.4f, %3.4f \n', uv_points(i, 1), uv_points(i, 2), uv_points(i, 3)); 
               fprintf('xyz : %3.4f, %3.4f, %3.4f \n', xyz_points(i, 1), xyz_points(i, 2), xyz_points(i, 3)); 
            end
                
            
            
        end
        
        function obj = addGentry(obj, xyz_world, uv_cam)
            
            u_cam = uv_cam(1) / uv_cam(3);
            v_cam = uv_cam(2) / uv_cam(3);
            
            Gentry = [xyz_world(1), xyz_world(2), xyz_world(3), 1, ...
                0, 0, 0, 0, ...
                -u_cam*xyz_world(1), -u_cam*xyz_world(2), -u_cam*xyz_world(3), -u_cam; ...
                0, 0, 0, 0, ...
                xyz_world(1), xyz_world(2), xyz_world(3), 1, ...
                -v_cam*xyz_world(1), -v_cam*xyz_world(2), -v_cam*xyz_world(3), -v_cam];
            
            
            if(size(obj.G_matrix,1) ==0)
                obj.G_matrix = Gentry;
            else
                obj.G_matrix = [obj.G_matrix; Gentry];
            end
        end
        
        function obj = constructPmatrix(obj)
            
            % create G matrix
            for i=1:size(obj.wld_points,1)
                obj= obj.addGentry(obj.wld_points(i,:), obj.cam_points(i, :));
            end
            
            % compute P matrix : the eigenvector of GtG with smallest
            % eigenvalue
            [U, S, V] = svd(obj.G_matrix);
            eigenv_gtg = V(:, end);
                    
            obj.proj_matrix_computed = [ eigenv_gtg(1) eigenv_gtg(2) eigenv_gtg(3) eigenv_gtg(4) ; ...
                            eigenv_gtg(5) eigenv_gtg(6) eigenv_gtg(7) eigenv_gtg(8) ; ...
                            eigenv_gtg(9) eigenv_gtg(10) eigenv_gtg(11) eigenv_gtg(12) ; ...
                            0                  0                  0                  1                  ];
            
        end
        
        function obj = getRotAndTfromPmatrix(obj)
            
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
            a

            r = inv(a) * B 
            [U, ~, V] = svd(r);
            d = sign(det(V*U'));
            r = [1 0 0; 0 1 0; 0 0 d]*r;

            t = inv(a) * b
            
            intr = c.intrinsic_camera_matrix(1:3, 1:3);
            t2 = inv(intr) * b
                     
            
            obj.R_noint = r;
            obj.t_noint = t;
            

        end
        
        
    end
    
    
end


classdef CamExtrinsicsCalibrationClass
    properties
        % initial camera matrix
        intrinsic_camera_matrix
        % extrinsics camera matrix (initial guess)
        extrinsic_camera_matrix
        % proj matrix initial guess
        proj_matrix
        
        cam_points
        cam_points_prior_noise
        wld_points
        
        % G matrix
        G_matrix
        %
        proj_matrix_computed
       
        R_ext
        t_ext
        
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
        
        function obj = CalculateExtrinsics(obj, wld_points, add_noise_bool)
            
            obj.wld_points = wld_points;
            
            obj.cam_points = (obj.proj_matrix*obj.wld_points')';
            obj.cam_points_prior_noise = obj.cam_points;
            
            snr = 10;
            if(add_noise_bool)
                 obj.cam_points = [awgn(obj.cam_points(:, 1:3), snr) obj.cam_points(:, 4)];
            end
                        
            obj = obj.constructPmatrix();
            
            obj = obj.getRotAndTfromPmatrix();
            
        end
        
        function obj = constructGmatrix(obj, cam_points, wld_points)
            
            if(isempty(cam_points) & isempty(wld_points))
                cam_points = obj.cam_points;
                wld_points = obj.wld_points;
            end
            
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
            
            obj.G_matrix = G_matrix;
        end
        
        function obj = constructPmatrix(obj)
            
            
                obj = obj.constructGmatrix([], []);
                % compute P matrix : the eigenvector of GtG with smallest
                % eigenvalue
                [U, S, V] = svd(obj.G_matrix);
                eigenv_gtg = V(:, end);

                pp = [ eigenv_gtg(1) eigenv_gtg(2) eigenv_gtg(3) eigenv_gtg(4) ; ...
                                eigenv_gtg(5) eigenv_gtg(6) eigenv_gtg(7) eigenv_gtg(8) ; ...
                                eigenv_gtg(9) eigenv_gtg(10) eigenv_gtg(11) eigenv_gtg(12) ; ...
                                0                  0                  0                  1                  ];

                pp = pp./(sqrt(pp(3,1)^2+pp(3,2)^2+pp(3,3)^2)) ;
                
                obj.proj_matrix_computed = pp;
        end
        
        function obj = getRotAndTfromPmatrix(obj)
            
            pp = obj.proj_matrix_computed;

            B = pp(1:3, 1:3); 
            b = pp(1:3, 4); 

            %             k = B*B';
            %             
            %             %u0
            %             a(1,3) = k(1,3); 
            %             %v0
            %             a(2,3) = k(2,3);
            %             %beta 
            %             a(2,2) = sqrt(k(2,2)-a(2,3)^2);
            %             %alpha
            %             a(1,1) = sqrt(k(1,1)-a(1,3)^2);
            % 
            %             r = inv(a) * B 
            %             [U, ~, V] = svd(r);
            %             d = sign(det(V*U'));
            %             r = [1 0 0; 0 1 0; 0 0 d]*r;
            % 
            %             t = inv(a) * b
            %             
            %             intr = obj.intrinsic_camera_matrix(1:3, 1:3);

            r = inv(obj.intrinsic_camera_matrix(1:3, 1:3)) * B;
            t = inv(obj.intrinsic_camera_matrix(1:3, 1:3)) * b;
            
            if(round(det(r)) == -1)
                [U, ~, V] = svd(r);
                d = sign(det(V*U'));
                r = -r;
                t = -t;
            end

            obj.R_ext = r;
            obj.t_ext = t;
            

        end
        
        
    end
    
    
end


clear
clc

path = '/home/cat/Documents/CMU_Herb/camera_ext_calibration/';
path = '~/cat_workspace/src/ExtrinsicsCalibration/matlab_code/';
addpath([path, 'ADAcode'])
addpath([path, 'bkg_code'])


nlevels = [0.02 0.05 0.1 0.3 ];
nr_experiments = 10;

% (1) Create world data from a table randomly 

n = 100;
tb_height = 0.5;

min_x = 0.5; max_x = 1;
min_y = -0.3; max_y = 0.3;
x = min_x + (max_x-min_x).*rand(n,1); 
y = min_y + (max_y-min_y).*rand(n,1); 
z = ones(n,1).* tb_height;
tag_xyz = [x y z ones(n, 1)];

for i=1:n
    tag_euler = [0 0 0.0001];
    tag_rotation = rotationAroundZ(tag_euler(3))* rotationAroundY(tag_euler(2))* rotationAroundX(tag_euler(1));
    T_w2tag(:, :, i) = [tag_rotation tag_xyz(i, 1:3)'; 0 0 0 1];
end

tag_size = 0.06;
ts = tag_size/2;
for i=1:n
    t = T_w2tag(:, :, i);
    c1 = t; c1(1,4) = t(1,4)+ts; c1(2,4) = t(2,4)+ts; 
    c2 = t; c2(1,4) = t(1,4)-ts; c2(2,4) = t(2,4)+ts;
    c3 = t; c3(1,4) = t(1,4)-ts; c3(2,4) = t(2,4)-ts;
    c4 = t; c4(1,4) = t(1,4)+ts; c4(2,4) = t(2,4)-ts;
    
    T_w2tag_c1(:, :, i) = c1; 
    T_w2tag_c2(:, :, i) = c2;
    T_w2tag_c3(:, :, i) = c3;
    T_w2tag_c4(:, :, i) = c4;
end

%% (2) Create a T_wc matrix
        min_ = -0.2; max_ = 0.2;
        rx = min_ + (max_-min_).*rand(1,1); 
        rz = min_ + (max_-min_).*rand(1,1); 

        R_init = rodrigues([rx 0 rz]);
        min_ = 0; max_= 1;
        t_init = [min_ + (max_-min_).*rand(1); min_ + (max_-min_).*rand(1); min_ + (max_-min_).*rand(1)];

        T_w2c_init = [R_init , t_init; 0 0 0 1];
        
for nl =1 : length(nlevels)
    for nexp=1 : nr_experiments
        
        
        noise_ratio = nlevels(nl);

        fprintf('Test nr %d\n', nexp);

        % using T_w2c (perfect matrix)
        tagposes = T_w2tag;
        for n =1:size(T_w2tag,3)
           % passing T_c_tag
           T_c_tag(:, :, n) = inv(T_w2c_init) * tagposes(:, :, n);
           T_c_c1(:, :, n) = inv(T_w2c_init) * T_w2tag_c1(:, :, n);
           T_c_c2(:, :, n) = inv(T_w2c_init) * T_w2tag_c2(:, :, n);
           T_c_c3(:, :, n) = inv(T_w2c_init) * T_w2tag_c3(:, :, n);
           T_c_c4(:, :, n) = inv(T_w2c_init) * T_w2tag_c4(:, :, n);
        end

        A = [529.2945         0  466.9604         0; ...
                    0  531.2835  273.2594         0; ...
                    0         0    1.0000         0; ...
                    0         0         0    1.0000];

        c1 = []; c2 = []; c3 = []; c4 = [];
        for n =1:size(T_w2tag,3)
           tc_1 = T_c_c1(:, :, n) ;
           corner1 = A * tc_1(1:4, 4);
           c1(n, :) = [corner1(1)/corner1(3) corner1(2)/corner1(3) 1];

           tc_2 = T_c_c2(:, :, n) ;
           corner2 = A * tc_2(1:4, 4);
           c2(n, :) = [corner2(1)/corner2(3) corner2(2)/corner2(3) 1];

           tc_3 = T_c_c3(:, :, n) ;
           corner3 = A * tc_3(1:4, 4);
           c3(n, :) = [corner3(1)/corner3(3) corner3(2)/corner3(3) 1];

           tc_4 = T_c_c4(:, :, n) ;
           corner4 = A * tc_4(1:4, 4);
           c4(n, :) = [corner4(1)/corner4(3) corner4(2)/corner4(3) 1];
        end


        % (3) Add Gaussian noise (statistically independent and addictive noise) to
        % data

        a = rodrigues(T_w2c_init(1:3, 1:3));
        N = randn(size(a));
        N = N*noise_ratio;
        a = a + N;
        R = rodrigues(a);
        t = T_w2c_init(1:3, 4);N = randn(size(t));N = N*noise_ratio;
        t = t + N;
        T_w2c_init_n = [R t; 0 0 0 1];

        counter = 0;
        for i=1:size(T_c_tag,3) 

            t_c = T_c_tag(:, :, i);
            t_w = T_w2c_init_n * t_c;

            zhat_tag = t_w(1:3, 3);
            z_tag =  t_w(3, 4);
            counter = counter + 1;
            norm_z(counter) = norm(z_tag - 0.5);

            CosTheta = dot(zhat_tag, [0 0 1]')/(norm(zhat_tag)*norm([0 0 1]'));
            theta_error(counter) = acos(CosTheta);

        end

        wld_points_train = T_w2tag(:, :, 1:n/2);
        wld_points_test = T_w2tag(:, :, n/2+1:end);
        cam_points_train = T_c_tag(:, :, 1:n/2);
        cam_points_test = T_c_tag(:, :, n/2+1:end);
        corners_train.c1 = c1(1:n/2, :);
        corners_train.c2 = c2(1:n/2, :);
        corners_train.c3 = c3(1:n/2, :);
        corners_train.c4 = c4(1:n/2, :);
        corners_test.c1 = c1(n/2+1:end, :);
        corners_test.c2 = c2(n/2+1:end, :);
        corners_test.c3 = c3(n/2+1:end, :);
        corners_test.c4 = c4(n/2+1:end, :);

        %% (4)run PB function

        %  find x
        w=rodrigues(T_w2c_init_n(1:3, 1:3));
        wx=w(1);
        wy=w(2);
        wz=w(3);

        x0 =[wx, wy, wz, T_w2c_init_n(1,4), T_w2c_init_n(2,4), T_w2c_init_n(3,4)];

        % optimization using fmincon
        options = optimoptions('fmincon', 'Algorithm', 'interior-point');
        [x, fval, exittag] = fmincon(@simulation_minimizeExtrinsics,x0,[],[],[],[],[],[],[], options, cam_points_train, corners_train, x0);

        % get refined extrinsics
        [R_est, t_est] = xToRt(x(1:6));

        T_w2c_PB = [R_est, t_est ; 0 0 0 1];

        %% Run World Known Poses Based Optimization (with test data)

        counter = 0;
        for i=n/2+1:size(T_c_tag,3)
            counter = counter + 1;
            cam_pose(:, counter) = [T_c_tag(:, 4, i)];
            wld_pose(:, counter) = [T_w2tag(:, 4, i)];
        end

        w=rodrigues(T_w2c_PB(1:3, 1:3));
        wx=w(1);
        wy=w(2);
        wz=w(3);

        x0 =[T_w2c_PB(1,4), T_w2c_PB(2,4), T_w2c_PB(3,4)];

        % optimization using fmincon

        options = optimoptions('fmincon','Algorithm', 'interior-point');

        [x, fval, exittag] = fmincon(@simulation_wldpose,x0,[],[],[],[],[],[],[], options, wld_pose, cam_pose, T_w2c_PB);

        % get refined extrinsics
        t_est = x(1:3)';

        T_w2c_WK = [T_w2c_PB(1:3, 1:3), t_est ; 0 0 0 1];

        %% Results

        % matrices 

        K_init(:,:, nexp) =T_w2c_init;

        K_init_n(:,:, nexp) =T_w2c_init_n;

        K_init_PB(:,:, nexp) =T_w2c_PB;

        K_init_WK(:,:, nexp) =T_w2c_WK;

        T_w2c_init
        T_w2c_init_n
        T_w2c_PB
        T_w2c_WK

        % world poses

        wld_pose;

        t_w_init_n = T_w2c_init_n *cam_pose;

        t_w_PB = T_w2c_PB*cam_pose;

        t_w_WK = T_w2c_WK*cam_pose;


        % error of projection 

        [mean_abs_error_init_n(nexp,:), std_deviation_init_n(nexp,:),  mean_error_3D_init_n(nexp), std_deviation_3D_init_n(nexp)]  = translationErrorBetweenPointsInWorld(wld_pose, t_w_init_n);

        [mean_abs_error_pb(nexp,:), std_deviation_pb(nexp,:),  mean_error_3D_pb(nexp), std_deviation_3D_pb(nexp)]  = translationErrorBetweenPointsInWorld(wld_pose, t_w_PB);

        [mean_abs_error_wk(nexp,:), std_deviation_wk(nexp,:),  mean_error_3D_wk(nexp), std_deviation_3D_wk(nexp)]  = translationErrorBetweenPointsInWorld(wld_pose, t_w_WK);


    end

    test_mean_init(nl,:) = mean(mean_abs_error_init_n)
    test_std_init(nl,:) = mean(std_deviation_init_n)

    test_mean_PB(nl,:) = mean(mean_abs_error_pb)
    test_std_PB(nl,:) = mean(std_deviation_pb)

    test_mean_WK(nl,:) = mean(mean_abs_error_wk)
    test_std_WK(nl,:) = mean(std_deviation_wk)

    test_mean_init_3D(nl) = mean(mean_error_3D_init_n)
    test_std_init_3D(nl) = mean(std_deviation_3D_init_n)

    test_mean_PB_3D(nl) = mean(mean_error_3D_pb)
    test_std_PB_3D(nl) = mean(std_deviation_3D_pb)

    test_mean_WK_3D(nl) = mean(mean_error_3D_wk)
    test_std_WK_3D(nl) = mean(std_deviation_3D_wk)

end

%% Plots


figure; 

x = [1,2,3,4]
scatter(x, test_mean_init_3D, 'r', 'MarkerFaceColor', 'r') 
hold on;
scatter(x, test_mean_PB_3D, 'g', 'MarkerFaceColor', 'g') 
hold on;
scatter(x, test_mean_WK_3D, 'b', 'MarkerFaceColor', 'b') 
title('Calibration errors with increased noise level')
legend('non Opt.', 'PB', 'WK')
%% 

red1 = [141 12 20]./255
red2 = [205 17 61 ]./255
red3 = [236 113 129]./255

blue1 = [40, 112, 220]./255
blue2 = [88 156 215]./255
blue3 = [75 178 222]./255

green1 = [13 94 61]./255
green2 = [14 138 88]./255
green3 = [75 222 163]./255



figure;
x = [1,2,3,4]
plot(x, test_mean_init(:, 1), '-s' , 'MarkerEdgeColor', red1 ,'MarkerFaceColor', red1) 
hold on;
plot(x, test_mean_PB(:, 1), '-s',  'MarkerEdgeColor', blue1,'MarkerFaceColor', blue1) 
hold on;
plot(x, test_mean_WK(:, 1), '-s',  'MarkerEdgeColor', green1,'MarkerFaceColor', green1) 
hold on;
plot(x, test_mean_init(:, 2), '-o', 'MarkerEdgeColor', red2 ,'MarkerFaceColor', red2 ) 
hold on;
plot(x, test_mean_PB(:, 2), '-o',  'MarkerEdgeColor', blue2,'MarkerFaceColor', 'c') 
hold on;
plot(x, test_mean_WK(:, 2), '-o',  'MarkerEdgeColor', green2,'MarkerFaceColor', green2) 
hold on;
plot(x, test_mean_init(:, 3), '-p',  'MarkerEdgeColor', red3 ,'MarkerFaceColor', red3) 
hold on;
plot(x, test_mean_PB(:, 3), '-p',  'MarkerEdgeColor', blue3,'MarkerFaceColor', blue3) 
hold on;
plot(x, test_mean_WK(:, 3), '-p',  'MarkerEdgeColor', green3,'MarkerFaceColor', green3) 

legend('non Opt., x', 'non Opt., y', 'non Opt., z','PB, x', 'PB, y', 'PB, z', 'WK, x', 'WK, y', 'WK, z')
clear
clc

% (1) Create world data from a table randomly 

n = 50;
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

% (2) Create a T_wc matrix 

T_w2c = [0.9754 0.2206 0 -0.0470; ...
         -0.2206 0.9754 0 0.05; ...
         0 0 1 -0.09; ...
         0 0 0 1];

% using T_w2c (perfect matrix)
tagposes = T_w2tag;
for n =1:size(T_w2tag,3)
   % passing T_c_tag
   T_c_tag(:, :, n) = inv(T_w2c) * tagposes(:, :, n);
   T_c_c1(:, :, n) = inv(T_w2c) * T_w2tag_c1(:, :, n);
   T_c_c2(:, :, n) = inv(T_w2c) * T_w2tag_c2(:, :, n);
   T_c_c3(:, :, n) = inv(T_w2c) * T_w2tag_c3(:, :, n);
   T_c_c4(:, :, n) = inv(T_w2c) * T_w2tag_c4(:, :, n);
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

ratio = 0.05;
a = rodrigues(T_w2c(1:3, 1:3));
N = randn(size(a));
N = N*ratio;
a = a + N;
R = rodrigues(a);
t = T_w2c(1:3, 4);N = randn(size(t));N = N*ratio;
t = t + N;
T_w2c_n = [R t; 0 0 0 1];

counter = 0;
for i=1:size(T_c_tag,3) 
        
    t_c = T_c_tag(:, :, i)
    t_w = T_w2c_n * t_c
  
    zhat_tag = t_w(1:3, 3);
    z_tag =  t_w(3, 4);
    counter = counter + 1;
    norm_z(counter) = norm(z_tag - 0.5);
    
    CosTheta = dot(zhat_tag, [0 0 1]')/(norm(zhat_tag)*norm([0 0 1]'));
    theta_error(counter) = acos(CosTheta);
       
end


%% (4)run function

%  find x

K_wc = T_w2c_n;

w=rodrigues(K_wc(1:3, 1:3))
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, K_wc(1,4), K_wc(2,4), K_wc(3,4)]

corners.c1 = c1;
corners.c2 = c2;
corners.c3 = c3;
corners.c4 = c4;
% optimization using fmincon
options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');
[x, fval, exittag] = fmincon(@simulation_minimizeExtrinsics,x0,[],[],[],[],[],[],[], options, T_c_tag, corners, x0);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1];

T_w2c

KK

T_w2c_n


%% error evaluation

for i=1:size(T_c_tag,3) 
        
    t_c = T_c_tag(:, :, i);
    t_w_estimated = KK * t_c;
    t_w_real = T_w2c * t_c;
    error_xyz(i) = sqrt(sum( (t_w_estimated(1:3, 4) - t_w_real(1:3, 4)).^2));
    
    error_rot(i) = sqrt(sum( (rodrigues(t_w_estimated(1:3, 1:3)) - rodrigues(t_w_real(1:3, 1:3)) ).^2));
    
end
error_xyz
sum(error_xyz)

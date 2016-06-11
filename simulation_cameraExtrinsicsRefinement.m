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

% (2) Create a T_wc matrix 

T_w2c = [0.9754 0.2206 0 -0.0470; ...
         -0.2206 0.9754 0 0.05; ...
         0 0 1 -0.09; ...
         0 0 0 1];


% (3) Add Gaussian noise (statistically independent and addictive noise) to
% data

ratio = 0.1;
a = rodrigues(T_w2c(1:3, 1:3));
N = randn(size(a));
N = N*ratio;
a = a + N;
R = rodrigues(a);
t = T_w2c(1:3, 4);N = randn(size(t));N = N*ratio;
t = t + N;
T_w2c_n = [R t; 0 0 0 1];

% for n =1:size(T_w2tag,3)
%     tt = T_w2tag(:, :, n);
%     
%     a = rodrigues(tt(1:3, 1:3));
%     N = randn(size(a));
%     N = N*ratio;
%     a = a + N;
%     %R = rodrigues(a);
%     R = tt(1:3, 1:3);
%     t = tt(1:3, 4);N = randn(size(t));N = N*ratio;
%     t = t + N;
%     T_w2tag_n(:, :, n) = [R t; 0 0 0 1];
% end

% using T_w2c (perfect matrix)
tagposes = T_w2tag;
for n =1:size(T_w2tag,3)
    T_c_tag(:, :, n) = inv(T_w2c) * tagposes(:, :, n);
end

counter = 0;
for i=1:size(T_c_tag,3) 
        
    t_c = T_c_tag(:, :, i)
    t_w = T_w2c_n * t_c
    
%     axis_angle = vrrotmat2vec(t_w(1:3, 1:3));
%     if(axis_angle(3) <0)
%         zhat_tag = -1 .* axis_angle(1:3);
%     else
%         zhat_tag = axis_angle(1:3);
%     end
%     
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

% optimization using fmincon

options = optimoptions('fmincon', 'Display','iter', 'Algorithm', 'interior-point');

[x, fval, exittag] = fmincon(@simulation_minimizeExtrinsics,x0,[],[],[],[],[],[],[], options, T_c_tag, x0);

% get refined extrinsics
[R_est, t_est] = xToRt(x(1:6));

KK = [R_est, t_est ; 0 0 0 1];

T_w2c

KK

T_w2c_n
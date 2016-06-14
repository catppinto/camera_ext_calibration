function [T] = performCalibration(w_points, c_points)


w_nomean = w_points - repmat(mean(w_points), size(w_points,1), 1) ; 
c_nomean = c_points - repmat(mean(c_points), size(c_points,1), 1);
% scatter3(w_nomean(:, 1), w_nomean(:, 2), w_nomean(:, 3), 'r');
% hold on;
% scatter3(c_nomean(:, 1), c_nomean(:, 2), c_nomean(:, 3), 'g');

%%
H = zeros(3,3);
for i=1:8
    H = H + ( w_nomean(i, :)'* c_nomean(i, :));
end

%%
[U, S, V] = svd(H); 

vut = V*U';
d = sign(det(vut)); 
R = V * [1 0 0 ; 0 1 0; 0 0 d] * U'

t = - R * mean(c_points)' + mean(w_points)'

T = [R t; 0 0 0 1];
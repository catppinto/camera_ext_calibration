function [w_p] = getWorldPointFromCameraPoint(c_p, T)
% Function discription : c_p must be a N*3 or N*4 point set. T is the
% transformation matrix to be used. w_p is the point in the world frame, a
% N*4 point / pointset


% R = T(1:3, 1:3);
% ww = zeros(8,3); cc = []; ww_icp =[];
% for i = 1:8
%     ww(i, :) = R * c_nomean(i, :)';
% end

if (size(c_p, 2) ~= 4) 
    c_p = [c_p ones(8,1)]';
end

ww = T * c_p;
w_p = ww';
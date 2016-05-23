function [error_, reproj_gt, reproj_test] = computeReprojectionError( points, gt, test)
% istance btw Non normalized uv points
n = size(points, 1); 

if(size(points, 2) ~= 4) 
    points= [points ones(size(points, 1), 1)];
end

reproj_gt = (gt * points')' ; 

reproj_test = (test * points')' ; 

% reproj_gt = [ reproj_gt(:, 1)./reproj_gt(:,3) reproj_gt(:, 2)./reproj_gt(:,3)];
% reproj_test = [ reproj_test(:, 1)./reproj_test(:,3) reproj_test(:, 2)./reproj_test(:,3)];
error_ =0;
for i=1:size(points,1)
    this_error = sqrt((reproj_gt(i,1)-reproj_test(i,1)).^2 + ...
        (reproj_gt(i,2)-reproj_test(i,2)).^2 + (reproj_gt(i,3)-reproj_test(i,3)).^2);
    error_ = error_ + this_error;
end
error_ = error_/size(points,1);
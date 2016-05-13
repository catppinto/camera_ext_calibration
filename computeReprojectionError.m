function error_ = computeReprojectionError( points, gt, test)

n = size(points, 1); 

if(size(points, 2) ~= 4) 
    points= [points ones(size(points, 1), 1)];
end

reproj_gt = (gt * points')' ; 

reproj_test = (test * points')' ; 

% reproj will be n X 4

error_ = sum(sum((reproj_gt(:, 1:3) - reproj_test(:, 1:3)).^2));
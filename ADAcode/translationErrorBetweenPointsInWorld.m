function [mean_abs_error, std_deviation,  mean_error_3D, std_deviation_3D] = translationErrorBetweenPointsInWorld(p1, p2)

n = size(p1, 2); 

error_xyz = p1 - p2;

mean_error = 1/n * sum(error_xyz,2);

mean_abs_error = 1/n * sum(abs(error_xyz),2);

std_deviation = sqrt(1/(n-1) * sum( error_xyz.^2 ,2));



% % the points should be in vector format !!! x; y; z
% error_t = sqrt((p1 - p2).^2);
% error_t = sum(error_t,2)/size(p1,2);
% %%

error = sqrt(sum((p1-p2).^2));
mean_error_3D = 1/n*sum(error, 2);

dev = error-mean_error_3D;
std_deviation_3D = sqrt(1/(n-1) * sum( dev.^2 ,2));



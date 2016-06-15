function [error_t, overall_error] = translationErrorBetweenPointsInWorld(p1, p2)
% the points should be in vector format !!! x; y; z
error_t = sqrt((p1 - p2).^2);
error_t = sum(error_t,2)/size(p1,2);
%%
overall_error = sqrt(sum( (p1 - p2).^2 ));
overall_error = sum(overall_error)/size(p1,2);

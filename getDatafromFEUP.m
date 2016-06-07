%% tag calibration data from FEUP : file markerarray_log

function [cam_points, wld_points] = getDatafromFEUP()

w_t0 = [0.7005,  0.024, 0]
w_t1 = [   0.4,  0.024, 0]
w_t2 = [0.2495,  0.024, 0]
w_t3 = [ 0.402,   0.19, 0]
w_t4 = [ 0.449, -0.004, 0.16627]
w_t5 = [ 0.449, -0.004, 0.11877]
w_t6 = [ 0.449, -0.004, 0.07127]
w_t7 = [ 0.449, -0.004, 0.02377]

uv_t0 = [65.4581, 254.323, 1]
uv_t1 = [365.963, 266.746, 1]
uv_t2 = [512.924, 273.093, 1]
uv_t3 = [375.592, 410.053, 1]
uv_t4 = [320.853, 93.2974, 1]
uv_t5 = [ 318.17, 142.573, 1]
uv_t6 = [315.855, 186.125, 1]
uv_t7 = [313.897, 225.133, 1]

cam_points = [uv_t0; uv_t1; uv_t2; uv_t3; uv_t4; uv_t5; uv_t6; uv_t7];

wld_points = [ w_t0;  w_t1;  w_t2;  w_t3;  w_t4;  w_t5;  w_t6;  w_t7];
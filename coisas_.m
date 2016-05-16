i = c.intrinsic_camera_matrix;
r = c.R_noint;

cp = c.cam_points(1,:);
wp = c.wld_points(1,:);

At = i(1:3, 1:3)*r*wp' - cp'


tz = At(3)/i(3,3)
tx = (At(1) - i(1,3)*tz)/a(1,1)
ty = (At(2) - i(2,3)*tz)/a(2,2)

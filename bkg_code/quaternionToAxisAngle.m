function axis_angle = quaternionToAxisAngle(q)
 
 w=q(1);
 x=q(2);
 y=q(3);
 z=q(4);
if (w > 1) 
   sqrtq= sqrt(w*w+x*x+y*y+z*z);
   q= q./sqrtq;
end

angle = 2 * acos(w);
s = sqrt(1-w*w); %assuming quaternion normalised then w is less than 1, so term always positive.
if (s < 0.001)  % test to avoid divide by zero, s is always positive due to sqrt
 % if s close to zero then direction of axis not important
 x = x; % if it is important that axis is normalised then replace with x=1; y=z=0;
 y = y;
 z = z;
else
 x = x / s; 
 y = y / s;
 z = z / s;
end

axis_angle = [angle, x, y, z];
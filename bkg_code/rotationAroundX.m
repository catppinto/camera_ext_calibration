function rot = rotationAroundX(s)

R_x =[1.0, 	 0.0, 	  0.0;
         0.0, cos(s), -sin(s);
         0.0, sin(s),  cos(s)];

rot = R_x;
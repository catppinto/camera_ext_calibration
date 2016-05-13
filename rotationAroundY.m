function rot = rotationAroundY(s)

    R_y =[ cos(s), 0.0, sin(s);
         0.0, 	 1.0, 	 0.0;
        -sin(s), 0.0, cos(s)];
rot = R_y;
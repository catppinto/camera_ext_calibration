function rot = rotationAroundZ(s)

 R_z =[ cos(s), -sin(s), 0.0;
         sin(s),  cos(s), 0.0;
            0.0,     0.0, 1.0];
rot = R_z;

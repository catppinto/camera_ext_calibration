function [R, t] = xToRt(x)

theta2=x(1)^2+x(2)^2+x(3)^2;
theta=sqrt(theta2);
omega=[0 -x(3) x(2);
    x(3) 0  -x(1);
    -x(2) x(1) 0;];
R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/(theta^2))*(omega*omega);
t =[x(4);x(5);x(6)];
function [R, t] = xToRt(x)

R = rodrigues(x(1:3));
t =[x(4);x(5);x(6)];
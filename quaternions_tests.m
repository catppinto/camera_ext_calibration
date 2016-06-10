%% Quaternions 


theta_deg = 0:1:360; 
theta = theta_deg./360 * (2*pi);
n =[0 0 1]'; 

% rotation around z axis of theta degrees
w = cos(theta/2);
x = n(1)*sin(theta/2);
y = n(2)*sin(theta/2);
z = n(3)*sin(theta/2);

q = [w' x' y' z']



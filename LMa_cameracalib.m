%%% https://engineering.purdue.edu/kak/computervision/ECE661/HW5_LM_handout.pdf



%% LMa for Camera Calibration Parameters Refinement

% au = alpha, bu = beta, sk = gamma
syms u0 v0 au av sk real
syms tx ty tz wx wy wz real
syms X Y um vm real

% intrinsics
K = [ au sk u0; ...
    0 av v0; ...
    0 0 1];


% expression for rotation matrix from rodrigues formula
theta = sqrt(wx^2 + wy^2 + wz^2);
omega = [ 0 -wz wy; ...
    wz 0 -wx; ...
    -wy wx 0];
R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);

% translation vector
t = [tx; ty; tz];

% perspective projection of the model point (X,Y)
uvs=K*[R(:,1) R(:,2) t]*[X; Y; 1];
u=uvs(1)/uvs(3);
v=uvs(2)/uvs(3);

% calculate the geometric distance in x and y direction
% um,vm =the x and y positions of an extracted corner
% u,v = the x and y positions of the projection of the corresponding model point
dx=um-u;
dy=vm-v;

% Evaluate the symbolic expression of the Jacobian w.r.t. the estimated parameters
Jx=jacobian(dx,[au,av,u0,v0,sk,wx wy wz tx ty tz]);
Jy=jacobian(dy,[au,av,u0,v0,sk,wx wy wz tx ty tz]);



%%

% R0, t0 = initial extrinsic parameters obtained from the solution in Section 3
tx=t0(1);
ty=t0(2);
tz=t0(3);

% convert the 3x3 rotation matrix into 3-vector w=[wx wy wz] of the Rodigrues representation
R=R0;
theta=acos((trace(R)-1)/2);
w=(theta/(2*sin(theta)))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
wx=w(1);
wy=w(2);
wz=w(3);


iters=200; 				% set the number of iterations for the LM algorithm
lamda=0.01; 				% initial the value of damping factor
updateJ=1;
Ndata=2*”number of corner points”
Nparams=11;
for it=1:iters
    if updateJ==1

        % create the intrinsic parameter matrix
        K=[au sk u0;
            0  av v0;
            0  0   1];

        % convert the 3-vector [wx wy wz] of the Rodigrues representation
        % into the 3x3 rotation matrix
        theta2=wx^2+wy^2+wz^2;
        theta=sqrt(theta2);
        omega=[0 -wz wy;
            wz 0  -wx;
            -wy wx 0;];
        R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
        t=[tx;ty;tz];

        % Evaluate the Jacobian at the current parameter values
        % and the values of geometric distance
        J=zeros(Ndata,Nparams);
        d=zeros(Ndata,1);
        for i=1:size(xe,2)
            X=Xw(1,i); 		% (X,Y) are the coordinates of the ith model point
            Y=Xw(2,i);
            J(2*(i-1)+1,:)= "the value of Jx evaluated at the ith model point X,Y and the current parameters"
            J(2*(i-1)+2,:)= "the value of Jy evaluated at the ith model point X,Y and the current parameters"

            % perspective projection of the model point
            uvs=K*[R(:,1) R(:,2) t]*[X; Y; 1];
            up=uvs(1)/uvs(3);
            vp=uvs(2)/uvs(3);
            % compute the geometric distance in x & y directions
            % xe(1,i), xe(2,i) = the the x and y positions of the ith extracted corner.
            %  up,vp = the x and y positions of the projection of the corresponding model point
            d(2*(i-1)+1,1)=xe(1,i)-up;
            d(2*(i-1)+2,1)=xe(2,i)-vp;
        end

        % compute the approximated Hessian matrix
        H=J'*J;
        if it==1 % the first iteration : compute the initial total error
            e=dot(d,d);
            disp(e);
        end
    end

    % Apply the damping factor to the Hessian matrix
    H_lm=H+(lamda*eye(Nparams,Nparams));

    % Compute the updated parameters
    dp=-inv(H_lm)*(J'*d(:));
    au_lm=au+dp(1);
    av_lm=av+dp(2);
    u0_lm=u0+dp(3)
    v0_lm=v0+dp(4)
    sk_lm=sk+dp(5);
    wx_lm=wx+dp(6);
    wy_lm=wy+dp(7);
    wz_lm=wz+dp(8);
    tx_lm=tx+dp(9);
    ty_lm=ty+dp(10);
    tz_lm=tz+dp(11);

    % Evaluate the total geometric distance at the updated parameters
    K=[au_lm sk_lm u0_lm;
        0    av_lm v0_lm;
        0    0       1];
    omega=[0 -wz_lm wy_lm;
        wz_lm 0  -wx_lm;
        -wy_lm wx_lm 0;];
    theta2=wx_lm^2+wy_lm^2+wz_lm^2;
    theta=sqrt(theta2);
    R = eye(3) + (sin(theta)/theta)*omega + ((1-cos(theta))/theta^2)*(omega*omega);
    t=[tx_lm;ty_lm;tz_lm];
    d_lm=zeros(Ndata,1);
    for i=1:size(xe,2)
        X=Xw(1,i);
        Y=Xw(2,i);
        uvs=K*[R(:,1) R(:,2) t]*[X; Y; 1];
        up=uvs(1)/uvs(3);
        vp=uvs(2)/uvs(3);
        d_lm(2*(i-1)+1,1)=xe(1,i)-up;
        d_lm(2*(i-1)+1,1)=xe(2,i)-vp;
    end

    % e_lm is the error between the image coordinates and projective coordinates using
    % the updated parameters
    e_lm=dot(d_lm,d_lm);

    %If the total geometric distance of the updated parameters is less than the previous one
    % then makes the updated parameters to be the current parameters
    % and decreases the value of the damping factor
    if e_lm<e
        lamda=lamda/10;
        au=au_lm;
        av=av_lm;
        u0=u0_lm;
        v0=v0_lm;
        sk=sk_lm;
        wx=wx_lm;
        wy=wy_lm;
        wz=wz_lm;
        tx=tx_lm;
        ty=ty_lm;
        tz=tz_lm;
        e=e_lm;
        disp(e);
        updateJ=1;
    else
        % otherwise increase the value of the damping factor and try again
        updateJ=0;
        lamda=lamda*10;
    end
end

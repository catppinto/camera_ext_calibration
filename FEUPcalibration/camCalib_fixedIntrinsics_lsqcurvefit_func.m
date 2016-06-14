function [res_lsq_noJacobian, res_lsq_Jacobian] = camCalib_fixedIntrinsics_lsqcurvefit_func(c)
% initial guess is now the result from Zhang paper

%% set initial guesses

t0 = c.t_ext;
R0 = c.R_ext;

%% Using matlab function lsqcurvefit

xdata = c.wld_points';
uv = c.cam_points';
ydata(1, :) = uv(1, :)./uv(3, :);
ydata(2, :) = uv(2, :)./uv(3, :);
 
% R0, t0 = initial extrinsic parameters obtained from the solution in Section 3
tx=t0(1);
ty=t0(2);
tz=t0(3);

w = rodrigues(R0);
wx=w(1);
wy=w(2);
wz=w(3);

x0 =[wx, wy, wz, tx, ty, tz];

display('xdata')
xdata
display('ydata') 
ydata
display('x0')
x0

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
lb = [];
ub = [];
[x1,resnorm,residual,exitflag,output]  = lsqcurvefit(@myfun,x0,xdata,ydata, lb,ub,options);
[R_est1, t_est1] = xToRt(x1);

res_lsq_noJacobian.x0 = x0;
res_lsq_noJacobian.ext0 = [R0 t0; 0 0 0 1];
res_lsq_noJacobian.options = options;
res_lsq_noJacobian.x = x1;
res_lsq_noJacobian.resnorm = resnorm;
res_lsq_noJacobian.residual = residual;
res_lsq_noJacobian.exitflag = exitflag;
res_lsq_noJacobian.output = output;
res_lsq_noJacobian.ext_est = [R_est1 t_est1; 0 0 0 1];

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'Jacobian','on');
lb = [];
ub = [];
[x2,resnorm2,residual2,exitflag2,output2] = lsqcurvefit(@myfun,x0,xdata,ydata, lb,ub,options);
[R_est2, t_est2] = xToRt(x2);

res_lsq_Jacobian.x0 = x0;
res_lsq_noJacobian.ext0 = [R0 t0; 0 0 0 1];
res_lsq_Jacobian.options = options;
res_lsq_Jacobian.x = x2;
res_lsq_Jacobian.resnorm = resnorm2;
res_lsq_Jacobian.residual = residual2;
res_lsq_Jacobian.exitflag = exitflag2;
res_lsq_Jacobian.output = output2;
res_lsq_Jacobian.ext_est = [R_est2 t_est2; 0 0 0 1];



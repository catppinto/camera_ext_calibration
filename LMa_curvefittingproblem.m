%% LM implementation on a line fitting problem 

syms a b x y real; 

f = (a * cos(b*x) + b * sin(a*x));

d = (y-f);
Jsym = jacobian(d, [a,b]);

%% Synthetic data generation with some noise
 
a = 100;
b = 102;
x = [0:0.1:2*pi]';
y = a*cos(b*x)+b*sin(a*x);
y_input = y + 5*rand(length(x), 1);

%% initial guess for parameters 
a0 = 100.5;
b0 = 102.5;

y_init = a0*cos(b0*x)+b0*sin(a0*x);
Ndata = length(y_input);
Nparams = 2;
n_iter = 100;
lamda = 0.01;
updateJ = 1;
a_est = a0 ; 
b_est = b0 ;
for it = 1:n_iter
    if updateJ == 1 
        %evaluate Jacobian matrix at current parameters 
        J = zeros(Ndata, Nparams);
        for i = 1:length(x)
            J(i,:) = [ - cos(b_est*x(i)) - b_est*x(i)*cos(a_est*x(i)), ...
                a_est*x(i)*sin(b_est*x(i)) - sin(a_est*x(i))];
        end
        % evaluate distance errors at current parameters 
        y_est = a_est*cos(b_est*x)+b_est*sin(a_est*x);
        d = y_input - y_est;
        %compute approximated Hessian 
        H = J' * J;
        if it == 1 %first iteration , compute total error
            e = dot(d,d);
        end
    end
    
    % applying damping factor on hessian
    H_lm = H + (lamda * eye(Nparams, Nparams));
    
    % Compute the updated parameters 
    dp=-inv(H_lm)*(J'*d(:)); 
    a_lm=a_est+dp(1); 
    b_lm=b_est+dp(2); 
    
    % Evaluate the total distance error at the updated parameters
    y_est_lm = a_lm * cos(b_lm*x) + b_lm * sin(a_lm *x); 
    d_lm=y_input-y_est_lm;     
    e_lm=dot(d_lm,d_lm); 
    % If the total distance error of the updated parameters is less than the previous one 
   % then makes the updated parameters to be the current parameters 
   % and decreases the value of the damping factor   
    if e_lm<e 
        lamda=lamda/10; 
        a_est=a_lm;        
        b_est=b_lm; 
        e=e_lm; 
        disp(e); 
        updateJ=1; 
    else % otherwise increases the value of the damping factor
        updateJ=0; 
        lamda=lamda*10; 
    end    
end








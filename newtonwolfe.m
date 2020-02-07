function x = newtonwolfe(p, x, i)


% Author      : Yuqing Chen
% Description : Newton's method with Wolfe line search
% Input       : p ~ problem function handle
%               x ~ initial iterate
%               i ~ input paramter value structure
% Output      : x ~ final iterate


% Set default value of i

% Store output strings
out_line = '============================';
out_data = '  k        F(x)        ||d||';
out_null =                   '----------';

% Print output header
fprintf('%s\n%s\n%s\n',out_line,out_data,out_line);

% Evaluate F at x
F = feval(p,x,0);

% Calculate the norm of the initial descent direction
g0 = max(1,norm(feval(p,x,1)));

% Initialize iteration counter
k = 0;


% Iteration loop
while 1
  
  % Print iterate information
  fprintf('%4d  %.4e  ',k,F);
  
  % Compute the descent firection
  g = -feval(p,x,1);  

  % Check termination conditions
  if k > i.maxiter || norm(g) <= i.opttol*g0, break; end
  
  % Evaluate Jacobian of F at x
  J = feval(p,x,1);
  
  % Evaluate Hessian of F at x
  H = feval(p,x,2);
  
  % Set initial value of modifing paramter
  xi = 1e-4;
  
  % Calculate the length of x
  n = length(x);
  
  % Modify H if H is not positive definite
  while min(eig(H))<0
      H = H + xi*eye(n);
      xi = 10*xi;
  end
  
  % Solve Newton system
  d = -H\J;
  
  % Evaluate norm of direction
  norms.d = norm(d);
  
  % Print search direction information
  fprintf('%.4e\n',norms.d);

  % Choose stepsize
  alpha = wolfe(p,x,i);
  
  % Update iterate
  x = x + alpha*d;
  
  % Evaluate F at x
  F = feval(p,x,0);
  
  % Evaluate error
  norms.F = norm(F);
  
  % Increment counter
  k = k + 1;
  
end

% Print final iterate information
fprintf('%s\n%s\n',out_null,out_line);

function x = steepestbacktrack(p, x, i)


% Author      : Yuqing Chen
% Description : Steepest descent with backtracking line search
% Input       : p ~ problem function handle
%               x ~ initial iterate
%               i ~ input paramter value structure
% Output      : x ~ final iterate



% Store output strings
out_line = '============================';
out_data = '  k        F(x)        ||d||';
out_null =                   '----------';

% Print output header
fprintf('%s\n%s\n%s\n',out_line,out_data,out_line);

% Initialize iteration counter
k = 0;

% Calculate initial function value
F = feval(p,x,0);

% Calculate initial descent direction
d = -feval(p,x,1);

% Calculate the norm of the initial descent direction
g0 = max(1,norm(d));
fprintf('%4d',g0);

% Iteration loop
while 1
  
  % Print iterate information
  fprintf('%4d  %.4e  ',k,F);
  
  % Compute the descent firection
  g = -feval(p,x,1);  
  
  % Check termination conditions
  if k > i.maxiter || norm(g) <= i.opttol*g0, break; end
  
  % Choose stepsize
  alpha = backtrack(p,x);

  % Update iterate
  x = x + alpha*g;
  
  % Evaluate functuion value
  F = feval(p,x,0);
  
  % Evaluate norm of direction
  norms.d = norm(feval(p,x,1));
  
  % Print search direction information
  fprintf('%.4e\n',norms.d);

  % Increment counter
  k = k + 1;
end

% Print final iterate information
fprintf('%s\n%s\n',out_null,out_line);
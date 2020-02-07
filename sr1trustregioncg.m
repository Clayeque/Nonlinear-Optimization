 function x = sr1trustregioncg(p, x, i)

% Author      : Yuqing Chen
% Description : SR1 quasi-Newton trust region method with CG subproblem solver
% Input       : p ~ problem function handle
%               x ~ initial iterate
%               i ~ input paramter value structure
% Output      : x ~ final iterate

% Set bound for trustregion method
Delta = 1;

% Evaluate F at x
F = feval(p,x,0);

% Calculate the norm of the initial descent direction
g0 = max(1,norm(feval(p,x,1)));

% Store output strings
out_line = '============================';
out_data = '  k        F(x)        ||d||';
out_null =                   '----------';

% Print output header
fprintf('%s\n%s\n%s\n',out_line,out_data,out_line);

% Initialize iteration counter
k = 0;

  
% Initialize Hessian of F at x
  H = feval(p,x,2);

while 1
  
  % Print iterate information
  fprintf('%4d  %.4e  ',k,F);
  
  % Compute the descent firection
  g = -feval(p,x,1);  
  

  % Check termination conditions
  if k > i.maxiter || norm(g) <= i.opttol*g0, break; end
  
  
  
  % Calculate the direction of Cauchy point
  dC = g'*g*g/(g'*H*g);
  
  % Calculate the direction of Newton point
  dN = cg(H,g,x,i);

  % Initinialize alpha
  a = 2;
  
  % Calculate alpha in Dogleg method
  if norm(dC) <= Delta
      a = 2;
  else
      if norm(dN) >= Delta
          a = Delta/norm(dN);
      else
      a = sqrt(Delta^2- norm(dC))/(norm(dC - dN))+1;
      end
  end
  
  if a >= 0 && a <= 1
      s = a*dC;
  else
      if a>=1 && a <=2 
        s = dC+(a-1)*(dN-dC);
      end
  end
  
  
  
  % Calculate actual reduction
  ar = feval(p,x,0)-feval(p,x+s,0);
  
  % Calculate predicted reduction
  pr = -feval(p,x,1)'*s - .5*s'*feval(p,x,2)*s;
  
  % Evaluate the direction
  rho = ar/pr;

  if rho >0
      x= x+s;
  end
  if rho > 0.75
      if norm(s)>0.8*Delta
          Delta = 2*Delta;
      end
  elseif rho<0.1
      Delta = 0.5*Delta;
  end
  
  
  % Update y
  y = feval(p,x,1) - g;

  if abs(s'*(y-H*s))>=i.sr1updatetol*norm(s)*norm(y-H*s)
      H = H + (s-H*y)*(s-H*y)'/((s-H*y)'*y);
  end
  
  
  % Evaluate norm of direction
  norms.s = norm(s);
  
  % Print search direction information
  fprintf('%.4e\n',norms.s);
  
  % Evaluate F at x
  F = feval(p,x,0);  
  
  % Evaluate error
  norms.F = norm(F);
  
  % Increment counter
  k = k + 1;
end


% Print final iterate information
fprintf('%s\n%s\n',out_null,out_line);


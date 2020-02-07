 function x = trustregioncg(p, x, i)

% Author      : Yuqing Chen
% Description : trust region method with CG subproblem solver
% Input       : p ~ problem function handle
%               x ~ initial iterate
%               i ~ input paramter value structure
% Output      : x ~ final iterate



% Set bound for trustregion method
Delta = 1;
Delta_hat = 1000;

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

while 1
  
  % Print iterate information
  fprintf('%4d  %.4e  ',k,F);
  
  % Compute the descent firection
  g = -feval(p,x,1);  
  
  % Check termination conditions
  if k > i.maxiter || norm(g) <= i.opttol*g0, break; end
  
  % Evaluate Hessian of F at x
  H = feval(p,x,2);
  
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
      d = a*dC;
  else
      if a>=1 && a <=2 
        d = dC+(a-1)*(dN-dC);
      end
  end
  
  
  % Calculate actual reduction
  ar = feval(p,x,0)-feval(p,x+d,0);
  
  % Calculate predicted reduction
  pr = -feval(p,x,1)'*d - .5*d'*feval(p,x,2)*d;
  
  % Evaluate the direction
  rho = ar/pr;

  if rho >0
        if rho < i.c1tr 
            Delta = Delta/4;
        else
            if rho>i.c2tr && abs(norm(d) - Delta)<1e-3
            Delta = min(2*Delta, Delta_hat);
            end
        end
        if rho > 1e-6
            x = x + d;
        end
  end
  % Evaluate norm of direction
  norms.d = norm(d);
  
  % Print search direction information
  fprintf('%.4e\n',norms.d);
  
  % Evaluate F at x
  F = feval(p,x,0);  
  
  % Evaluate error
  norms.F = norm(F);
  
  % Increment counter
  k = k + 1;
end

% Print final iterate information
fprintf('%s\n%s\n',out_null,out_line);



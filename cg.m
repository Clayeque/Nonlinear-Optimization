function x = cg(A,b,x,i)

% function x = cg(A,b,x,e)
%
% Author      : Frank E. Curtis
% Description : Conjugate Gradient method (preliminary version).
% Input       : A ~ symmetric matrix
%               b ~ right-hand side vector
%               x ~ initial iterate
%               cgopttol ~ solution tolerance
%               cgmaxiter ~ iteration limit
% Output      : x ~ final iterate
% Revised by  : Yuqing Chen



% Evaluate initial iterate norm
norms.x = norm(x);

% Evaluate initial residual
r0 = A*x-b;

% Evaluate residual norm
norms.r = norm(r0);

% Store initial residual norm
norms.r0 = norms.r;

% Evaluate initial direction
p = -r0;

% Evaluate direction norm
norms.p = norm(p);

% Initialize iteration counter
k = 0;

% Main CG loop
while 1
  
 
  % Check for termination
  if norms.r <= i.cgopttol*max(norms.r0,1) || k > i.cgmaxiter, break; end
  
  % Evaluate matrix-vector product
  Ap = A*p;
  
  % Evaluate vector-vector product
  pAp = p'*Ap;
  
  % Evaluate steplength
  alpha = (r0'*r0)/pAp;
  
  % Update iterate
  x = x + alpha*p;
  
  % Evaluate residual
  r1 = r0 + alpha*Ap;
  
  % Evaluate residual norm
  norms.r = norm(r1);
  
  % Evaluate CG multiplier
  beta = (r1'*r1)/(r0'*r0);
  
  % Update direction
  p = -r1 + beta*p;
  
  % Evaluate direction norm
  norms.p = norm(p);
  
  % Renew residual
  r0 = r1;
  
  % Increment iteration counter
  k = k + 1; 
  
end

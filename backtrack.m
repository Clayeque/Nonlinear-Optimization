function alpha = backtrack(p,x)

% Compute stepsize with backtracking line search
% Input       : p ~ problem function handle
%               x ~ initial iterate
% Output      : alpha ~ stepsize with backtracking line search

  
  % Set the initial value of alpha
  alpha = 1;
  
  % Compute the descent firection
  g = -feval(p,x,1); 
  
  while 1
      
      % Check Armijo-Goldstein conditions
      if feval(p, x+alpha*g, 0) <= feval(p,x,0) + 0.5*alpha*feval(p,x,1)'*g
          break;
      else
          alpha = alpha/2;
      end
  end

end
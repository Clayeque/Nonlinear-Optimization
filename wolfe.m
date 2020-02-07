function alpha = wolfe(p,x,i)

% Compute stepsize satisfying the Wolfe condition
% Input       : p ~ problem function handle
%               x ~ initial iterate
% Output      : alpha ~ stepsize satisfying the Wolfe condition
  

  % Set the upper bound of alpha
  alpha_max=1;
  
  % Set the initial value of alpha
  alpha = 1;
  
  % Compute the descent firection
  g = -feval(p,x,1); 
  
  while 1
      
      % Check Wolfe conditions
      if feval(p, x+alpha*g, 0) <= feval(p,x,0) + i.c1ls*alpha*feval(p,x,1)'*g
          if feval(p,x+alpha*g,1)'*g >= i.c2ls*feval(p,x,1)'*g
          break;
          else
              alpha = (alpha+alpha_max)/2;
          end
      else
           alpha = alpha/2;  
      end
  end
end
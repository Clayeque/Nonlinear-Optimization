function [ x ] = cgTrust(p, x0, toler,initdel,maxdel,eta)
%  Implements the Steihaug-Toint conjugate gradient trust region method for
%  finding an approximate solution to the subproblem:
%
%    min m(p) = f + g'p + 1/2 * p' B p        s.t. ||p|| <= Del
%  Input:
%    fun      - a pointer to a function
%    x0       - starting point         
%  Output:
%    x      - the solution structure, with the solution point along with
%             function, gradient, and Hessian evaluations there of.
% Set gradient tolerance.
% Set initial delta value.
% Set maximum delta value.
% Set eta value.
del = initdel;  % Set delta value to initial delta value.
maxit = 1000;  % Set maximum number of allowed iterations.
% Evaluate initial iterate norm
xk=x0;
norms.x0 = norm(x0);
Ak=feval(p,xk,2);
bk=feval(p,xk,1);
% Evaluate initial residual
rk = Ak*xk-bk;
% Evaluate residual norm
norms.rk = norm(rk);
% Store initial residual norm
norms.r0 = norms.rk;
% Evaluate initial direction
p0 = -rk;
% Evaluate direction norm
norms.p0 = norm(p0);
% Initialize iteration counter
k = 0;
for k = 0:maxit
%  Compute  gradient at current point.
    gk = feval(p, xk, 1);
    %  Check for termination condition: norm of gradient less than toler.
    if norm(gk) < toler
        x = xk;
    end   
% Conjugate gradient method.
    while 1
        % Check for termination
        if norms.r <= e*max(norms.r0,1) 
            break; 
        end
        % Evaluate matrix-vector product
        Apk = A*pk;
        % Evaluate vector-vector product
        pkApk = pk'*Apk;
        % Evaluate steplength
        alpha = -(rk'*pk)/pkApk;
        % Update iterate
        xk = xk + alpha*pk;
        % Evaluate iterate norm
        norms.xk = norm(xk);
        % Evaluate residual
        rk = Ak*xk - bk;
        % Evaluate residual norm
        norms.rk = norm(rk);
        % Evaluate CG multiplier
        beta = (rk'*Apk)/pkApk;
        % Update direction
        pk = -rk + beta*pk;
        % Evaluate direction norm
        norms.pk = norm(pk);
        % Increment iteration counter
        k = k + 1;
    end
    %  Compute the reduction ratio.
    rho = (feval(p,xk,0) - feval(p,xk+pk,0))/(-gk.g'*pk - 0.5*pk.'*Ak.*pk);
    %  Update the trust region; i.e., update del and the current point.
    if rho < 0.25
        del = del / 4;
    else
        %  Note that radius only increases if norm of p reaches the trust
        %  region boundary.
        if rho > 0.75 && norm(p) == del
            del = min(2*del, maxdel);
        end
    end
    if rho > eta
        xk = xk + p;
    end
end
x=xk;
end

% function p = boundary(p, q, del)
% a = norm(q)^2;
% b = 2*p'*q;
% c = norm(p)^2 - del^2;
% alpha = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
% p = p + alpha * q;
% end
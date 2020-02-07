function x=nlpsolver(p,x,a,i)
% Author      : Yuqing Chen
% Description : Nonlinear optimization programming solver with different
% algorithms
% Input       : p ~ problem function handleis
%               x ~ initial iterate
%               i ~ input paramter value structure
% Output      : x ~ final iterate
tic;
% Set deta paramter structure
if ~isfield(i,'opttol')
    i.opttol = 1e-6;
end

if ~isfield(i,'maxiter')
    i.maxiter = 1e+3;
end
if ~isfield(i,'c1ls')
    i.c1ls = 0.1;
end
if ~isfield(i,'c2ls')
    i.c2ls = 0.9;
end
if ~isfield(i,'c1tr')
    i.c1tr = 0.25;
end
if ~isfield(i,'c2tr')
    i.c2tr = 0.75;
end
if ~isfield(i,'cgopttol')
    i.cgopttol = 1e-6 ;
end
if ~isfield(i,'cgmaxiter')
    i.cgmaxiter = 1e+3;
end

if ~isfield(i,'sr1updatetol')
    i.sr1updatetol = 1e-8;
end

if ~isfield(i,'bfgsdamptol')
    i.bfgsdamptol = 0.2;
end



if strcmp(a, 'steepestbacktrack')==1
    x=steepestbacktrack(p,x,i);
    
elseif strcmp(a, 'steepestwolfe')==1
    x=steepestwolfe(p,x,i);
    
elseif strcmp(a, 'newtonbacktrack')==1
    x=newtonbacktrack(p,x,i);
    
elseif strcmp(a, 'newtonwolfe')==1
    x=newtonwolfe(p,x,i);
    
elseif strcmp(a, 'trustregioncg')==1
    x=trustregioncg(p,x,i);
    
elseif strcmp(a, 'sr1trustregioncg')==1
    x=sr1trustregioncg(p,x,i);
    
elseif strcmp(a, 'bfgsbacktrack')==1
    x=bfgsbacktrack(p,x,i); 
    
elseif strcmp(a, 'bfgsbackwolfe')==1
    x=bfgswolfe(p,x,i); 

end
toc;
end

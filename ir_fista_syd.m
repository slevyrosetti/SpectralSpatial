function [xs, info] = ir_fista_syd(x, varargin)
%function [xs, info] = ir_fista_syd(x, varargin)
%|
%| Generic FISTA method with cost calculation
%| Adapted by Sydney Williams, Michigan 2015
%|
%| in
%|	x	[np 1]		initial estimate
%|
%| required
%|	'gradfun' @()		function returning gradient: gradfun(x)
%|	'step'			user-selected step size (1/Lipschitz)
%|
%| option
%|	'proxfun' @()		function performing proximal step
%|					(default: @(x)=x)
%|					e.g. @(x) = max(x,0)
%|					for nonnegativity constraint
%|  'costfun' @()       function calculating cost at each iteration
%|                  (default: @(x)=x)
%|                  to reduce runtime
%|  'frac'              integer fraction to set cost stopping criterion
%|                  (default 10)
%|  'norm_tol'          stopping tolerance for change in relative norm
%|                  ||x(n+1)-x(n)||/||x(0)|| (default 1e-3)
%|  'max_tol'           stopping tolerance for change in max abs value
%|                  ||x(n+1)-x(n)||_inf/||x(0))|| (default 1e-5)
%|  'design_tol'        stopping tolerance for relative norm to target pattern
%|                    ||x(n+1)-x(n)||/||d.*W|| (default 1e-4)
%|  'stop'      0|1|2|3	set to 0 to eliminate stopping criterion (default: 1)
%|                  1: norm change 2: max value change 3: norm design change
%|	'niter'			# total iterations (default 1)
%|	'chat'			verbosity (default 0)
%|	'restart'	0|1	1 to use restart (default: 1)
%|	'momentum'	0|1	set to 0 to eliminate momentum (default: 1)
%|  'target' [nz 1] design target pattern
%|  'weight' [nz 1] mask used to weight target pattern
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter 6]	iteration, cost, time, l-2 norm change, l-inf norm change, l-2 norm change to target pattern
%|
%| Copyright 2015-07-01, Jeff Fessler, University of Michigan
%| Updated   2015-07-28, 2015-08-10, 2015-10-16, 11-29-16 Sydney Williams, University of Michigan

if nargin == 1 && streq(x, 'test')
	run_mfile_local('ir_fista_test'), return
end
if nargin < 3, help(mfilename), error(mfilename), end

arg.gradfun = [];
arg.proxfun = @(x) x; 
arg.costfun = @(x) x;                 % filler, not really the cost
arg.frac = 10;
arg.norm_tol=1e-3;
arg.max_tol=1e-5;                     
arg.stop = 1;               
arg.niter = 1;
arg.step = [];
arg.chat = 0;
arg.isave = [];
arg.restart = true;
arg.alpha_restart = -cos(80*pi/180); % 80 degrees
arg.momentum = 1;
arg.weight = 1;
arg.target = 1;
arg.design_tol=1e-4;

arg = vararg_pair(arg, varargin);
arg.isave = iter_saver(arg.isave, arg.niter);


if isempty(arg.gradfun), fail('gradfun required'), end
if isempty(arg.step), fail('step required'), end

if ~isa(arg.gradfun, 'function_handle'), error 'gradfun not function handle?', end
if ~isa(arg.proxfun, 'function_handle'), error 'proxfun not function handle?', end

xs = zeros(numel(x), length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

info = zeros(arg.niter, 6);

% fast iterative shrinkage/thresholding algorithm (FISTA)
v = x;
norm_int=norm(x);                       % Added S. Williams 08-10-15
max_int=norm(x,inf);
xold = x; told = 1;
cost_max=arg.costfun(x);
cpu etic
for iter=1:arg.niter
	ticker(mfilename, iter, arg.niter)

	grad = arg.gradfun(v);

	x = v - arg.step * grad;
	x = arg.proxfun(x);

	% adaptive restart: todo - why not working?
	tmp = ((v - x)' * (x - xold)) / norm(v - x) / norm(x - xold);
%	pr tmp
	if arg.restart && (tmp > arg.alpha_restart)
		t = 1;
		v = x;
		if arg.chat
			printm('restart at %d', iter)
		end
	else
		t = (1 + sqrt(1 + 4 * told^2)) / 2;
		frac = arg.momentum * (told-1) / t;
		v = x + frac * (x - xold);
    end
    cost_old=arg.costfun(xold);             % Added S. Williams 08-10-15
    norm_change=norm(xold-x)/norm_int;      % Added S. Williams 08-10-15
    max_change=norm(xold-x,inf)/...
        max_int;                            % Added S. Williams 10-16-15
    norm_design=norm(xold-x)/...            % Added S. Williams 11-29-16
        norm(arg.target.*arg.weight);
	xold = x; told = t;
    
	if any(arg.isave == iter)               % Calculates cost at each iteration
		xs(:, arg.isave == iter) = x;
        cost(iter)=arg.costfun(x);          % Added S. Williams 07-28-15
    end
    
    % prints out info: 1-iteration number 2-cost at iteration 3-cpu time
    % 4-change in norm difference 5-change in max difference
    info(iter,1) = iter;                    % Added S. Williams 07-28-15
    info(iter,2) = cost(iter);              % Added S. Williams 07-28-15
    info(iter,3) = cpu('etoc');             % Added S. Williams 07-28-15
    info(iter,4) = norm_change;             % Added S. Williams 08-10-15   
    info(iter,5) = max_change;              % Added S. Williams 10-16-15
    info(iter,6) = norm_design;             % Added S. Williams 11-29-16

    % imposes stopping criterion: 1-norm change difference 2-max difference
    if arg.stop~=0
        if arg.stop==1 && norm_change<arg.norm_tol  % Stopping criterion
            xs=xs(:,1:iter);                        % Added S. Williams 08-10-15
            info=info(1:iter,:);
            return
        elseif arg.stop==2 && max_change<arg.max_tol% Added S. Williams 10-16-15
            xs=xs(:,1:iter);
            info=info(1:iter,:);
            return
        elseif arg.stop==3 && norm_design<arg.design_tol%Added S. Williams 11-29-16
            xs=xs(:,1:iter);    
            info=info(1:iter,:);
            return
        else
        end
    end
          
end
	
end

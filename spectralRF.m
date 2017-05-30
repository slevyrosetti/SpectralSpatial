function [b1,gx,gy,gz,m_approx,wn,Mxy,A,W] = spectralRF(d,fmap,roi,Trf,lambda,zslice,f)
%function [b1,gx,gy,gz, m_approx, wn, Mxy, A, W] = spectralRF(d,fmap,roi,Trf,lambda, zslice,f)
% Design tip-down or tip-up pulse using spectrally selective pulse
% Modified from Hao Sun, University of Michigan 2013 by Sydney Williams 2016

% INPUTS:
% d                [nx ny nz]   - prewinding spectral target pattern
% fmap             [nx ny nz]   - 3D fieldmap
% roi              [nx ny nz]   - 3D logical mask
% Trf              [1]          - pulse length in s
% lambda           [1]          - regularization term, typically set to 0
% zslice           [1]          - slice number used for 2D design
% f                [nf|0]       - if non-zero, empirical design target bandwidth used for
%                                 spectral prewinding; if zero, use full off-resonance bandwidth
% OUTPUTS:
% b1               [npts 1]    - spectral RF pulse in G
% gx               [npts 1]    - excitation gradient x (zeros)
% gy               [npts 1]    - excitation gradient y (zeros)
% gz               [npts 1]    - excitation gradient z (zeros)
% m_approx         [nf 1]      - small-tip angle approximation of magnetization in rad
% wn               [nf 1]      - vector of frequencies used in design in Hz
% Mxy              [nf 1]      - spectral prewinding target pattern
% A                [nf npts]   - MRI system matrix for small-tip angle pulse
% W                [nf nf]     - diagnoal weighting matrix (identity)

if Trf > 20e-3
   error('Trf is in s');
end
%% Construct frequency response
TEeff = median(angle(d(roi))/2/pi./(fmap(roi)+0.001));  % !!! there maybe phase wrap problem here
mxy_abs = median(abs(d(roi)));
fmap_roi = fmap(:,:,zslice);
fmap_roi = fmap_roi(roi(:,:,zslice)); 
if f==0                                                 % use full bandwidth range
    wn = min(fmap_roi):max(fmap_roi);
    wn = wn'; % convert to column vector
else
    f_av=mean(fmap_roi);                                % use emprical bandwidth range
    wn=f+f_av;
    wn=wn';
end
Mxy = mxy_abs*exp(1i*2*pi*wn*TEeff);                    % spectral prewinding target pattern

%% Construct System Matrix A
m0 = 1;                                                 % initial magnetization vector
gambar = 4257;                                          % gyromagnetic ratio in G/Hz
gam = gambar*2*pi;                                      
dt = 4e-6;                                             	% dwell time in s; 
t = 0:dt:Trf;
t = t-Trf;                                              % time-reversed vector for prewinding
npts=length(t);                                         % number of time points in RF pulse
[b0tt, tt] = ndgrid(wn(:),t(:));        
A = 1i*gam*m0*exp(1i*2*pi*(b0tt.*tt))*dt;               % spectral system matrix
    
%% Design RF pulse b1
cmax=0.2;                                               % Peak amplitude constraint in Gauss
cnorm=0.75;                                             % Norm power constraint  
ntol=1e-6;                                              % Norm change stopping criterion
niter_cg=150;                                           % # iterations used for intialization w/conj. grad.
niter_gp=1000;                                          % max # iterations used with FISTA
ntol=1e-6;                                              % Norm change stopping criterion
mtol=(2^-12)*(cmax)*1e-4;                               % Max value change stopping criterion, modified to be in terms of T
asplit=0.5;                                             % Split factor for proximal operator, if used

    
b1=zeros(npts,1);                                       % place holdder for RF pulse 
% Iterative CG pulse
[b1, info] = qpwls_pcg1(b1, A, ...                      % use unreguarlized CG solution for initialization of FISTA
    (eye(length(wn),length(wn))), Mxy(:),...
    sqrt(lambda), 'niter', niter_cg);
sys=A'*A;                                               % Need to find largest eigenvalue for step size of gradient projection
[vec,val]=powerit(ones(npts,1),sys,1e-5);               % Power iteration method to find largest eigenvalue
stepsize=1/val;
func.grad = @(x) -(A'*(Mxy(:)-A*x));                    % Negative gradient
func.peak = @(x) min(abs(x),cmax).*exp(sqrt(-1)*...
angle(x));                                              % Peak amplitude constraint, maintains design phase
func.power = @(x) ((min(cnorm-norm(x),0)<0)*...
(cnorm*(x/norm(x)))+(min(cnorm-norm(x),0)>=0)*x);       % Power constraint, related to SAR though not exact
func.peak_power= @(x) (asplit)*(min(abs(x),cmax).*...
exp(sqrt(-1)*angle(x)))+(1-asplit)*(x/(max(...
(norm(x)/cnorm),1)));                                   % Proximal average of peak amplitude and power constraints
func.cost = @(x) (0.5)*(Mxy(:)-A*x)'*(Mxy(:)-A*x);      % Cost function
% FISTA constrained pulse
[xs,info]=ir_fista_syd(b1, 'gradfun', func.grad,'costfun',func.cost, 'step', stepsize, 'restart', 0, ...
'proxfun',func.power,'stop',2,'norm_tol',ntol,'max_tol',mtol,'niter', niter_gp, 'isave', 'all');
b1=xs(:,end);                                           % final iteration
gx = zeros(size(b1));                                   % no excitation gradients
gy = zeros(size(b1)); 
gz = zeros(size(b1));
W=eye(length(wn),length(wn));                           % spectral "weights" (none for spectral)
kx=gx; ky=gy; kz=gz;                            
m_approx=A*b1;                                          % small-tip angle approximation of magnetization (in radian)

return
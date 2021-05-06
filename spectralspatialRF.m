function [b1, gx,gy,gz, m_approx, f, d, circ, kx,ky,kz, A, W, time] = spectralspatialRF(d,b0map,roi,Trf,FOV_x,FOV_y,underx,undery,zslice,targetf,ktype,varargin)
% function [b1, gx,gy,gz, m_approx, f, d, circ, kx,ky,kz, A, W] = spectralspatialRF(d,b0map,roi,Trf,FOV_x,FOV_y,underx,undery,zslice,targetf,ktype,varargin)
% Design spectral-spatial prewinding tip down pulse with 2D spatial inhomogenity correction using a pre-loaded fieldmap
% By Sydney Williams, University of Michigan, 2016
%
% INPUTS:
% (required)
% d                [nx ny nz]   - prewinding spectral target pattern
% b0map            [nx ny nz]   - 3D fieldmap in Hz
% roi              [nx ny nz]   - 3D logical mask
% Trf              [1]          - pulse length in ms
% FOV_x            [nx]         - FOV in x-direction in cm, usually same as FOV_y
% FOV_y            [ny]         - FOV in y-direction in cm, usually same as FOV_x
% underx           [1]          - spatial down-sampling factor in x direction
% undery           [1]          - spatial down-sampling factor in y direction
% zslice           [1]          - slice number used for 2D design
% targetf          [1]          - local target tracking bandwidth in kHz
% ktype            [1-9]        - excitation trajectory type, see excite_vds_grad.m for descriptions
% (optional)
% 'beta'           [1]          - regularization term, default = 0
% 'xoff'           [0|1]        - turn off x-gradients for 1D design, default = 0
% 'yoff'           [0|1]        - turn off y-gradients for 1D design, default = 0
%
% OUTPUTS:
% b1               [npts 1]    - spectral RF pulse in G
% gx               [npts 1]    - excitation gradient x (zeros)
% gy               [npts 1]    - excitation gradient y (zeros)
% gz               [npts 1]    - excitation gradient z (zeros)
% m_approx         [nf nx ny]  - small-tip angle approximation of magnetization in rad
% f                [nf 1]      - vector of frequencies used in design in Hz
% d                [nf nx ny]  - spectral-spatial prewinding target pattern
% A                [nf*nx*ny npts]- MRI system matrix for small-tip angle pulse
% W                [nf nx ny]  - diagnoal spectral-spatial weighting matrix 

fm=b0map;
arg.beta=0;                                         % default regularization, none
arg.xoff=0;                                         % default x gradients are on
arg.yoff=0;                                         % default y gradients are on
arg=vararg_pair(arg, varargin);                     % reads in user-provided values
TEeff = median((angle(d(roi))/2/pi)./(fm(roi)+0.001)); % !!! there maybe phase wrap problem here
mxy_abs = median(abs(d(roi)));


gambar = 42570;                                     % gamma/2pi in kHz/T
gam = gambar*2*pi;                                  % gamma in 1000*radians/s/T
dt = 4e-3;                                          % time sample size in ms (4 microseconds)
tfree = 2*TEeff*1e3;                                % free precession time in ms
tau = Trf*1e3;                                       % length of pulse in ms
time = [0:dt:tau]';                                 % time frame ms
npts = length(time);                                % number of time points for simulation

%% Excitation Trajectory Spatial Specifications
                 
gamp=5.0;                           % peak grad amp in G/cm 
nil=1;                              % number of interleaves
smax=15;                            % slew rate max in G/cm/ms (18 standard, 15 to reduce PNS)
% returns gradient, k-space, and slew waveforms for given trajectory type 'ktype'
[gx,gy,gz,kx,ky,kz,sx sy,sz] = excite_vds_grad(ktype,FOV_x,dt,gamp,smax,npts,nil);
if arg.xoff
    gx=zeros(npts,1);               % turns off gradients in x
    kx=zeros(npts,1);
end
if arg.yoff
    gy=zeros(npts,1);               % turns off gradients in y
    ky=zeros(npts,1);
end
% defines spatial sample space as well as downsamples ROI
[nobj_x, nobj_y, slices]=size(fm);                  % Spatial samples in x, y, and number of slices y of field map
nobj_x=ceil(nobj_x/underx);                         % # X samples, undersampled
nobj_y=ceil(nobj_y/undery);                         % # Y samples, undersampled
pixel_x=FOV_x/nobj_x;                               % pixel size in x (cm)
pixel_y=FOV_y/nobj_y;                               % pixel size in y (cm)
x=[-nobj_x/2:(nobj_x/2-1)]*pixel_x;                 % excitation space in x-direction
y=[-nobj_y/2:(nobj_y/2-1)]*pixel_y;                 % excitation space in y-direction
[X,Y]=ndgrid(x,y);                                  % vector for designing field map mask
bmap=fm(1:underx:end,1:undery:end,zslice)*1e-3;     % downsamples field map, converts to kHz
circ=roi(1:underx:end,1:undery:end,zslice);         % downsamples ROI

%% Spectral-spatial Design Setup
% df=10*1e-3;                                         % sampling rate, 0.01kHz or 10Hz
% wtype=1;                                            % type of weighting matrix (1 trad. specspat, 2 varying weight specspat, 3 purely spectral)
% shape=1;                                            % shape of weighting matrix used (1 rect.m 2 Gaussian spread)
% [W,f,b0map_tp,Nmap, b0map_val]=weight_matrix...
%     (bmap,circ,df,targetf,wtype,shape);             % creates weighting matrix and sampled frequency vector
% minf=min(f); maxf=max(f); nobj_f=length(f);         % descriptor variables about sampled frequencies
% [F, X, Y]=ndgrid(f, x, y);                          % grids frequency and spatial samples
%spsp=[F(:) X(:) Y(:)].';                            % stacks two row vectors   
t=time-tau;                                         % backwards time for spectral design
%kspace=[t kx ky];                                   % stacks spectral time and k-space trajectory into two column vectors
% d=mxy_abs*exp(1i*2*pi*f*tfree/2);                   % target design d
% d=repmat(d.',[1 nobj_x nobj_y]);                    % replicates spectral design matrix over length of spatial domain

% =========================================================================
% SLR: design 3D target pattern (fXxXy)for fat saturation pulse
nom_fa=90;                                          % required flip angle in degrees
f_fat=-3.3*127.74;                                  % fat frequency at 3T
minf=1.1*f_fat;                                     % take margins
maxf=100;                                           % include water frequency
df=10;                                              % sampling rate of 10Hz
f=minf:df:maxf;                                     % frequency axis
nobj_f=length(f);                                   % number of sampling frequencies
[~, ~, Y]=ndgrid(f, x, y);                          % grids frequency and spatial samples
exc_BW=50;                                          % excitation (for fat) and non-excitation (for water) BW of 50Hz
d=zeros(nobj_f, nobj_x, nobj_y);                   % placeholder for the 3D target pattern
d((f_fat-exc_BW/2) < f & f < (f_fat+exc_BW/2), :, :)=sind(nom_fa);  % 90 degree flip angle aroung fat and 0 elsewhere
% W=zeros(nobj_f, nobj_x, nobj_y);                   % placeholder for the weighting matrix
% W((f_fat-exc_BW/2) < f & f < (f_fat+exc_BW/2), :, :)=1; % optimize pulse around the fat frequency
% W((-exc_BW/2) < f & f < (exc_BW/2), :, :)=1;            % optimize pulse around the water frequency
W=ones(nobj_f, nobj_x, nobj_y);
% =========================================================================

N=size(Y);                                          % for use in Gnufft fatrix object
nshift_x=N(2)/2;                                    % centers X dimension
nshift_y=N(3)/2;                                    % centers Y dimension
deltaf=(maxf-minf)/(nobj_f-1); 
nshift_f=-minf/deltaf;                              % centers f dimension (not necessarily at 0 Hz) 
Nshift=[nshift_f nshift_x nshift_y];                % for NUFFT
om=2*pi*[t*deltaf, kx*pixel_x, ky*pixel_y];         % for use in Gnufft fatrix object
%minmax(om)
A=Gnufft({om, N, [6 6 6], 2*N, Nshift, 'table', 2^12, 'minmax:kb'});
A=1i*gam*dt*A';                                     % Conversion to proper units/dimensions (reverse from reconstruction)

%% Iterative Spectral-spatial RF Pulse Design

b1=zeros(npts,1);                                   % placeholder for RF pulse
cnorm=0.75e-04;                                     % max RF integrated power (T^2*ms)
cmax=0.2e-04;                                       % max peak RF power (T)
niter_cg=150;                                       % # iterations used for intialization w/conj. grad.
niter_gp=1000;                                      % max # iterations used with FISTA
ntol=1e-6;                                          % norm change stopping criterion
mtol=(2^-12)*(cmax);                                % max value change stopping criterion 
asplit=0.5;                                         % split for prox. operator
    
[b1, info] = qpwls_pcg1(b1, A, ...                  % creates initialization pulse    
    Gdiag(W), d(:), sqrt(arg.beta), 'niter', niter_cg);  
sys=A'*Gdiag(W)*A;                                  % need to find largest eigenvalue for step size of gradient projection
[vec,val]=powerit(ones(npts,1),sys,1e-5);           % power iteration method to find largest eigenvalue
stepsize=1/val;                                     % estimate of Lipschitz constant
func.grad = @(x) -(A'*(Gdiag(W)*(d(:)-A*x)));       % negative gradient
func.peak = @(x) min(abs(x),cmax).*exp(sqrt(-1)*...
angle(x));                                      	% Peak amplitude constraint, maintains design phase
func.power = @(x) ((min(cnorm-norm(x),0)<0)*...
(cnorm*(x/norm(x)))+(min(cnorm-norm(x),0)>=0)*x); 	% Power constraint, related to SAR though not exact
func.peak_power= @(x) (asplit)*(min(abs(x),cmax).*...
exp(sqrt(-1)*angle(x)))+(1-asplit)*(x/(max(...
(norm(x)/cnorm),1)));                               % Proximal average of peak amplitude and power constraints          
func.cost = @(x) (0.5)*(d(:)-A*x)'*Gdiag(W)...
    *(d(:)-A*x);                                    % Cost function
[xs,info]=ir_fista_syd(b1, 'gradfun', func.grad,'costfun',func.cost, 'step', stepsize, 'restart', 0, ...
'proxfun',func.power,'stop',2,'norm_tol',ntol,'max_tol',mtol,'niter', niter_gp, 'isave', 'all');
b1=xs(:,end);

m_approx=reshape(A*b1, [nobj_f nobj_x nobj_y]);     % small-tip angle approximation of magnetization (in radian)
b1=b1*1e4;                                          % converts from T to G

return; 

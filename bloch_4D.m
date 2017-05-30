function [mBloch, fBloch, b0map_tp] = bloch_4D(b1, gx, gy, gz, b0map, roi, simuRange, Tread, varargin)
% function [mBloch, fBloch] = bloch_4D(b1, gx, gy, gz, b0map, roi, simuRange, Tread, varargin)
% Computes 4D Bloch simulation over 3 spatial dimensions x,y,z as well as
% intravoxel frequency range f
% Sydney Williams, University of Michigan 2016
%
% INPUTS:
% (required)
% b1               [npts nc]     - RF pulses (potentially multi-transmit)
% gx               [npts 1]      - excitation x gradient
% gy               [npts 1]      - excitation y gradient
% gz               [npts 1]      - excitation z gradient
% b0map            [nx ny nz]    - field map
% roi              [nx ny nz]    - binary ROI mask
% simuRange: .x    [nx]          - contains x, y, z samples
%            .y    [ny]
%            .z    [nz]
% Tread            [1]           - free precession time of sequence in ms (2x TE)
%
% (optional)
% 'dt'             [1]           - gradient dwell time (ms), default = 4e-3
% 'sens'           [nc nx ny nz] - sensitivity maps (all ones if single transmit)
% 'zslice'         [1]           - z-th slice for 2D design
% 'wtype'          [1|2]         - type of bandwidth spread for simulation, default = 1 
%                                  (1: same across all x,y 2: spatially-varying w/off-resonance)
% 'shape'          [1|2]         - shape of weights for bandwidth spread, default = 1 
%                                  (1: rect spread, all 1's for all frequencies
%                                   2: gaussian spread of weights across f between 0 and 1 for +/- 3 st.dev.)
% 'bwmin'          [1]           - minimum intravoxel bandwidth spread for simulation (Hz), default = 25 Hz
% 'fsamp'          [1]           - sampling rate for intravoxel frequencies (Hz), default = 10 Hz
%
% OUTPUTS:
% mBloch           [nx ny nz nf] - Bloch simulation (3D if not including z-dimension)
% fBloch           [nf]          - frequency vector used in simulation
% b0map_tp         [nx ny nz]    - spatially varying intravoxel bandwidths (2D if not including z-dimension)

% sets optional arguments
[npts,nc]=size(b1);                 % number of time points and number of transmit coils
[nx,ny,nz]=size(b0map);             % number of x, y, z samples
arg.dt=4e-3;                        % default dwell time in ms
arg.sens=ones(nc,nx,ny,nz);         % default sensitivity maps
arg.wtype=1;                        % default bandwidth type for simulation
arg.shape=1;                        % default weights for frequencies in bandwidth
arg.bwmin=25;                       % default minimum bandwidth in Hz
arg.fsamp=10;                       % default frequency sampling rate in Hz
arg.zslice=[];                      % default 2D slice (none, 3D)
arg = vararg_pair(arg, varargin);   % reads in user-provided values

if isempty(arg.zslice)              % adjusts field map, roi, and sens. maps for 2D vs 3D
    b0mapz=b0map;
    roiz=roi;
    sens=arg.sens;
else
    b0mapz=b0map(:,:,arg.zslice);
    roiz=roi(:,:,arg.zslice);
    simuRange.z=0;
    sens=arg.sens(:,:,:,arg.zslice);
end

% sets up frequency range for simulation
[Wsim,fsim,b0map_tp,Nmap]...             	% determines spatially varying intravoxel bandwidth (b0map_tp) and 
    =weight_matrix(b0mapz,...            	% frequency varying weights across bandwidth Wsim
    roiz,arg.fsamp,arg.bwmin,arg.wtype,arg.shape);
gauss_sigma=max(max(b0map_tp));             % Gaussian st. dev. (max intravoxel spread)
minf=-3*gauss_sigma;
maxf=3*gauss_sigma;                         % simulate over +/- 3 st. dev.'s

mode=0;                                   	% simulate after one TR
T1=Inf; T2=Inf;                             % set to infinity for simulation
necho=round(Tread/2/arg.dt)                  	% number of time samples at echo time (halfway through free precession)
b1e = [b1; zeros(necho,nc)];                % simulate the excitation pattern at the echo time 
gxe = [gx; zeros(necho,1)];        
gye = [gy; zeros(necho,1)];
gze = [gz; zeros(necho,1)];
% 4D Bloch simulation
[mBloch, ~,fBloch] = parallel_Bloch4D(0,b1e,gxe,gye,gze,sens,simuRange.x,...
    simuRange.y,simuRange.z,arg.dt,fliplr(b0mapz),fliplr(roiz),T1,T2,minf,maxf,arg.fsamp,mode);
fprintf('done.');

return

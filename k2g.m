function g = k2g(k_kpcm, fov_ispcm, dt, gambar)
%function g = k2g(k_kpcm, fov_ispcm, dt, gambar)
% Converts k-space trajectory (cm^-1) to gradient trajectory (G/cm)
% Hao Sun, University of Michigan 2013
% Commented by Sydney Williams, University of Michigan 2016
% INPUTS:
% (required)
% k_kpcm           [npts nd]    - k-space trajectory trajectory (number of dimensions)
%    - k                          cycle/fov, require fov_ispcm be [Nd]
%    - kpcm cycle/cm              require strcmp(fov_ispcm, 'pcm') == true
% fov_ispcm        [nd]         - FOV relationship for k-space, 1 cm
%    - ispcm       'pcm'          string
% (optional)
% 'dt'             [1]          - dwell time in ms
% 'gambar'         [0|1]        - gyromagnetic ratio in 2pi*kHz/T
%
% OUTPUTS:
%  g               [Nt, Nd]     - gradient trajectory, G/cm 
% Gradient from kTraj, modified from Hao's k2g_hao().
% INPUT
%  k_kpcm [Nt, Nd]
%    - k,    cycle/fov, require fov_ispcm be [Nd]
%    - kpcm, cycle/cm,  require strcmp(fov_ispcm, 'pcm') == true
%  fov_ispcm
%    - fov,   [Nd]    , cm
%    - ispcm, 'pcm'   , string
%  dt,  [1]     , optional, mSec
%  gambar, [1]  , optional, 2pi*kHz/T
% OUTPUT
%  g,   [Nt, Nd], G/cm

if ~exist('dt',     'var'), dt = 4e-3; end % unit in mSec
if ~exist('gambar', 'var'), gambar = 2*pi*42.576e3; end % 2pi * kHz/T

dtSec = dt*1e-3;        % convert to Sec to pair with Hz
gamma = gambar/2/pi/10; % convert to Hz/Gauss

if strcmp(fov_ispcm, 'pcm') == true
  kpcm = k_kpcm;
elseif isnumeric(fov_ispcm)
  fov = fov_ispcm;
  k = k_kpcm;
  kpcm = bsxfun(@rdivide, [k; k(end,:)], fov);
end

g = diff(kpcm,1,1)/gamma/dtSec;
end

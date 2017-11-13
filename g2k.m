function [k, kpcm] = g2k(g, fov, dt, gambar)
%function [k, kpcm] = g2k(g, fov, dt, gambar)
% Converts gradient trajectory (G/cm) to k-space trajectory (cm^-1)
% Hao Sun, University of Michigan 2013
% Commented by Sydney Williams, University of Michigan 2016
% INPUTS:
% (required)
% g                [npts nd]    - gradient trajectory (number of dimensions 1-3), G/cm
% FOV              [nd]         - FOV relationship for k-space, 1 cm
% (optional)
% 'dt'             [1]          - dwell time in ms
% 'gambar'         [0|1]        - gyromagnetic ratio in 2pi*kHz/T
%
% OUTPUTS:
%  k               [Nt, Nd]     - k-space trajectory, FOV/cm 
%  kpcm            [Nt, Nd]     - k-space trajectory, cm^-1

if ~exist('dt',     'var'), dt = 4e-3; end % unit in mSec
if ~exist('gambar', 'var'), gambar = 2*pi*42.576e3; end % 2pi * kHz/T

dtSec = dt*1e-3;                % convert to Sec
gamma = gambar/2/pi/10;         % convert to Hz/Gauss

kpcm = -flipud(cumsum(flipud(g),1) * dtSec * gamma);
k = bsxfun(@times, kpcm, fov);

end
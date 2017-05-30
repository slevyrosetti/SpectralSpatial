function [gw]=gauss_weight(f,sigma)
% creates Gaussian distribution to be used for weighted average of frequency simulations for intravoxel dephasing
% Sydney Williams, University of Michigan 2016
% 
% Inputs:
% f         [nf]        vector of simulation frequencies (Hz) 
% sigma     [1]         Gaussian value of sigma (Hz) from Mulkern et al. 2016
%
% Outputs:
% gw        [nf]        normalized Gaussian weights for use in weighted average


    sigmarad=sigma*2*pi;    % Converts Gaussian width parameter to rad
    frad=f*2*pi;            % Converts simulation frequenices to rad

    gw=(1/sigmarad*(2*pi)^(0.5))*exp(-(frad).^2/(2*sigmarad^2));
    gw=gw/max(gw);          % normalizes to 1

end
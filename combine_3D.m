function [meanFBloch, abs_meanFBloch, angle_meanFBloch]=combine_3D(mBloch, f, b0map_tp, fsamp, zslice)
% function [meanFBloch, abs_meanFBloch, angle_meanFBloch]=combine_4D(mBloch, f, b0map_tp, fsamp)
% Combines 4D Bloch simulated magnetization into the z-th 2D slice with
% intravoxel frequenices combined used a (Gaussian) weighted distribution
%  Sydney Williams, University of Michigan 2016
%
% INPUTS:
% mBloch           [nx ny nz nf] - Bloch simulation for x, y, z spatial and f spectral domains 
%                                  (could just be x,y, and f)
% f                [nf]          - full frequency range used in 4D Bloch simulation
% b0map_tp         [nx ny nz]    - target bandwidth frequencies (Hz) at each spatial location, 
%                                  may account for through-plane gradient weights or be uniform
% samp_f           [1]           - frequency sampling rate in Hz
% zslice           [1]           - z-th slice number to evaluate

% OUTPUTS:
% meanFBloch       [nx ny]]      - complex average magnetization for 2D slice
% abs_meanFBloch   [nx ny]]      - average magnitude magnetization for 2D slice
% angle_meanFBloch [nx ny]]      - average phase magnetization for 2D slice

if ndims(mBloch)<3 || ndims(mBloch)>4, fail('must be 2 or 3 spatial dimensions'), end

if ndims(mBloch)==4                 % selects zth slice if 3 spatial dimensions
    mBloch=squeeze(mBloch(:,:,zslice,:));
end
[nx ny nf]=size(mBloch);            % size of multi-frequency simulation
gwmat=zeros(nx,ny,nf);              % place holder Gaussian matrix for weights
for i=1:nx
    for j=1:ny
        gauss_sigma=b0map_tp(i,j);  % gaussian standard deviation at each x,y location          
        fw=freq_vector(-3*gauss_sigma,3*gauss_sigma,fsamp); % frequencies for weights
        gw=gauss_weight(fw,gauss_sigma); % computes weights
        gw=gw/sum(gw);              % normalizes to sum to 1
        if isnan(gw)                % fills in NAN
            gw=0;
        else                        % finds starting location for Gaussian weights
            fstart=find(abs(f-...   % finds closest frequency index to start at
                min(fw))==min(abs(f-min(fw))));
            gwmat(i,j,fstart(1)...  % fills in weights to Gaussian matrix
                :(fstart(1)+length(gw)-1))=gw;
        end
    end
end
meanFBloch=sum(mBloch.*gwmat,3)...  % weighted mean complex magnetization
    ./sum(gwmat,3);
meanFBloch(isnan(meanFBloch))=0;   	% fills in NAN
abs_meanFBloch=abs(meanFBloch);
angle_meanFBloch=angle(meanFBloch)*180/pi;
return



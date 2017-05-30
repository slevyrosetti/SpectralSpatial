function [W,f,b0map_tp,Nmap,b0map_val] = weight_matrix(b0map,roi,df,targetf,type,shape)
% function [W,f]= function weight_matrix(b0map,roi,df,targetf,type)
% Create spectral-spatial weighting matrix based on range of frequencies
% from fieldmap
% Sydney Williams, University of Michigan, 2017
%
% INPUTS:
% b0map     [Nx Ny]     - 2D fieldmap in Hz
% roi       [Nx Ny] 	- 2D ROI associated with fieldmap
% df        [1]         - frequency sample rate in Hz
% targetf 	[1]         - minimum target bandwidth (intravoxel frequency spread) in Hz
% type      [1]         - 1: spectral-spatial case with same target bandwidth L
%                         everywhere (25 Hz used in first specspat paper submit)
%                         2: spectral-spatial case with varying target bandwidth L
%                         based on through-plane gradient from Yip 2009
%                         3: purely spectral pulse (no weights)
% shape     [1|2]       - 1: Include rect (all ones) weighting of samples included in design
%                       - 2: Include Gaussian weights for selected samples
%
% OUTPUTS:
% W         [Nf Nx Ny]  - 3D spectral-spatial weighting matrix
% f         [Nf]        - frequency samples vector
% b0map_tp  [Nx Ny]     - target bandwidth frequencies at each spatial location, 
%                         may account for through-plane gradient weights or be uniform
% Nmap      [Nx Ny]     - map of number of frequency samples included in
%                         weighting matrix at each spatial location
% b0map_val [Nx Ny]     - 2D fieldmap in Hz centered at the median value
%                        
gambar = 42570;                             % gamma/2pi in kHz/T
mask=b0map.*roi;
[sval,sind]=sort(b0map(roi),'descend');   	% sorts mask and will take 2nd largest and 2nd smallest freq in field map
[nobj_x,nobj_y]=size(b0map);                % spatial dimensions
mid_val=median(b0map(roi));                 % median value of the fieldmap (offset) in Hz   
b0map_val=(b0map(roi)-mid_val);             % values of fieldmap within ROI centered to 0 Hz
if abs(sval(1))<df && abs(sval(end))<df     % if values are within sampling frequency range,      
   maxf=sval(1)+targetf/2;                	% highest target frequency in kHz
   minf=sval(end)-targetf/2;              	% lowest  target frequency in kHz
   f=minf:maxf;                             % sample every 1Hz
   nobj_f=length(f);                        % length of frequency fector
   W=ones(nobj_f,nobj_x,nobj_y);            % weighting matrix includes all frequencies
else
   if type==1                              	% original method with fixed target bandwidth
        if shape==1
            maxf=sval(1)+targetf/2;         % highest target frequency + L/2 in kHz 
            minf=sval(end)-targetf/2;       % lowest  target frequency - L/2 in kHz
        elseif shape==2
            maxf=sval(1)+3*targetf;         % highest target frequency +3 sigma in kHz
            minf=sval(end)-3*targetf;       % lowest  target frequency -3 sigma in kHz
        end
        f=minf:df:maxf;                 	% increments frequency range by sampling rate from min to max value
        nobj_f=length(f);                   % length of frequency fector
        W=zeros(nobj_f,nobj_x,nobj_y);  	% placeholder for weighting matrix  
        if shape==1
        	nfreq=round(targetf/(f(2)-f(1)));	% number of frequency samples needed to cover bandwidth
        elseif shape==2
           	nfreq=round(6*targetf/(f(2)-f(1)));	% number of frequency samples needed to cover bandwidth
        end
        offres=zeros(nobj_x,nobj_y);      	% placeholder for off resonance values at each location
        b0map_tp=targetf*double(roi);       % local target bandwidth at each spatial location
        Nmap=nfreq*double(roi);             % number of frequency samples included in weighting matrix at each spatial location
        for i=1:nobj_x
            for j=1:nobj_y
                [~,center]=min(abs(f-mask(i,j)));  	% finds center resonance frequency at spatial location
                low=center-round(nfreq/2);         	% lower bound of excitation bandwidth
                high=center+round(nfreq/2);        	% upper bound of excitation
                offres(i,j)=center;                	% off-resonance value at each location
                if low<1                         	% corrects for wrapping around, elongates excitation bandwidth on alternate side
                    high=high+abs(low);
                    low=1;
                elseif high>nobj_f
                    low=low-(high-nobj_f);
                    high=nobj_f;
                end
                fvec=[(low:high)-median(low:high)]*df;
                gw=gauss_weight(fvec,targetf);  % Gaussian weights
                if roi(i,j)>0                   % adds weight to frequency bandwidth for spatial location
                    if shape==1
                        W(low:high,i,j)=1;    	% adds rect weight to frequency bandwidth for spatial location
                    elseif shape==2
                        W(low:high,i,j)=gw;     % adds Gaussian weight to frequency bandwidth for spatial location
                    end
                end
            end
        end
   elseif type==2
        dbdz=-2e-4;                         % throughplane gradient from Yip'09, G/cm/Hz
        zthick=0.4;                         % slice thickness in cm
        b0map_grad=(dbdz*zthick*...
            gambar/10)*b0map_val;           % b0map throughplane variation
        b0map_width=sqrt((b0map_grad).^2+...% local bandwidth L adjustment based on throughplane gradient
            targetf^2);
        b0map_tp=embed(b0map_width,roi);    % embeds bandwidths in 2D slice ROI for throughplane gradient
        minbw=min(b0map_width);             % minimum target bandwidth
        maxbw=max(b0map_width);             % maximum target bandwidth
        if shape==1
            maxf=sval(1)+maxbw/2;           % highest target frequency in kHz
            minf=sval(end)-maxbw/2;         % lowest  target frequency in kHz
        elseif shape==2
            maxf=sval(1)+3*maxbw;           % highest target frequency in kHz
            minf=sval(end)-3*maxbw;         % lowest  target frequency in kHz
        end
        f=minf:df:maxf;                 	% increments frequency range by sampling rate from min to max value
        nobj_f=length(f);                   % length of frequency fector
        W=zeros(nobj_f,nobj_x,nobj_y);  	% placeholder for weighting matrix 
        offres=zeros(nobj_x,nobj_y);      	% placeholder for off resonance values at each location
        for i=1:nobj_x
            for j=1:nobj_y
                if shape==1
                     nfreq(i,j)=round(b0map_tp(i,j)...
                     /(f(2)-f(1)));          % number of frequency samples needed to cover bandwidth
                elseif shape==2
                     nfreq(i,j)=round(6*b0map_tp(i,j)...
                     /(f(2)-f(1)));          % number of frequency samples needed to cover bandwidth
                end
                [~,center]=min(abs(f-mask(i,j)));% finds center resonance frequency at spatial location
                low=center-round(nfreq(i,j)/2);% lower bound of excitation bandwidth
                high=center+round(nfreq(i,j)/2);% upper bound of excitation
                offres(i,j)=center;        	% off-resonance value at each location
            	if low<1                   	% corrects for wrapping around, elongates excitation bandwidth on alternate side
                    high=high+abs(low);
                    low=1;
                elseif high>nobj_f
                    low=low-(high-nobj_f);
                    high=nobj_f;
                end
                fvec=[(low:high)-median(low:high)]*df;
                gw=gauss_weight(fvec,b0map_tp(i,j));
                if roi(i,j)>0            	% adds weight to frequency bandwidth for spatial location
                    if shape==1
                        W(low:high,i,j)=1;	% adds Gaussian weight to frequency bandwidth for spatial location
                    elseif shape==2
                        W(low:high,i,j)=gw; % adds rect weight to frequency bandwidth for spatial location
                    end
                end
            end
        end
        Nmap=nfreq;                         % number of frequency samples included in weighting matrix at each spatial location
   elseif type==3
        maxf=sval(1)+targetf/2;            	% highest target frequency in kHz
        minf=sval(end)-targetf/2;       	% lowest  target frequency in kHz
        f=minf:df:maxf;                   	% sample every 1Hz
        nobj_f=length(f);                	% length of frequency fector
        W=ones(nobj_f,nobj_x,nobj_y);     	% weighting matrix all ones (spectral pulse)
        b0map_tp=ones(nobj_x,nobj_y);       % no spectral target bandwidth
       	Nmap=ones(nobj_x,nobj_y);          
    end
end
b0map_val=embed(b0map_val,roi);             % embeds target bandwidth values within ROI
end
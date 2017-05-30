function [btd,gxtd,gytd,gztd,btu,gxtu,gytu,gztu,d,m_atEcho,m_beforeTipup,m_afterTipup] = prewinding_pulse(pulseType, seqType, b0map, roi, simuRange, zslice, nom_fa, Tread, varargin)
% function [btd,gxtd,gytd,gztd,btu,gxtu,gytu,gztu,d,m_atEcho,m_beforeTipup,m_afterTipup] = prewinding_pulse(pulseType, seqType, b0map, roi, simuRange, zslice, nom_fa, Tread, varargin)
% Designs prewinding RF pulses (purely spectral or spectral-spatial)
% Sydney Williams, University of Michigan 2016
%
% INPUTS:
% (required)
% pulseType        [1|2|3]      - prewinding pulse type, 1: purely spectral from Hao Sun 2016
%                                 2: 2D spectral-spatial from Sydney Williams, 2017
% seqType          [1|2]        - sequence type  1: STFR 2: SPGR
% b0map            [nx ny nz]   - field map
% roi              [nx ny nz]  	- binary ROI mask
% simuRange: .x    [nx]        	- contains x, y, z samples
%            .y    [ny]
%            .z    [nz]
% zslice           [1]          - z-th slice for 2D design
% nom_fa           [1]          - target flip angle (degrees)
% Tread            [1]           - free precession time of sequence in ms (2x TE)
%
% (optional)
% 'targetf'        [1]          - target local bandwidth for spectral-spatial pulse, default = 25 Hz
% 'ktype'          [0-10]       - excitation trajectory for spectral-spatial pulse, default = 9 (2 alternating VDS)
%                                 1: 1VDS at end 2: 1VDS at middle 3:1VDS at end 
%                                 4: 2VDS low-res even 5: 2VDS low-res even alternating
%                                 6: 3VDS even 7: 3VDS eve nalternating 8: 2VDS med-res
%                                 9: 2VDS med-res alternating 10: 1VDS high-res 0: purely spectral pulse
% 'undersamp'      [1]          - undersampling spatial factor, default = 4
% 'spec_band'      [0|1]        - spectral bandwidth for purely spectral pulse, default = 0
%                                 0: full bandwidth (Sun 2016) 1: bandwidth narrowed to empirical limit (Williams 2017)
% 'Trf'            [1]          - RF pulse length in sec, default = 3e-3
% 'sens'           [nc nx ny nz]- sensitivity maps (all ones if single transmit)
% 'doplot'         [0|1]        - 1: plot magnetization at TE and after tip-up (if STFR) 0: do not plot, default = 0
% OUTPUTS:
% btd              [npts]       - tip-down prewinding pulse
% gxtd             [npts]       - excitation x-gradient for tip-down pulse
% gytd             [npts]       - excitation y-gradient for tip-down pulse
% gztd             [npts]       - excitation z-gradient for tip-down pulse
% btu              [npts]       - tip-up prewinding pulse
% gxtu             [npts]       - excitation x-gradient for tip-up pulse
% gytu             [npts]       - excitation y-gradient for tip-up pulse
% gztu             [npts]       - excitation z-gradient for tip-up pulse
% d                [nx ny]      - 2D simulated magnetization at TE
% m_atEcho         [nx ny]      - 2D simulated magnetization at TE
% m_beforeTipup    [nx ny]      - 2D simulated magnetization before tip-up pulse (set to 0 for SPGR)
% m_afterTipup     [nx ny]      - 2D simulated magnetization after tip-up pulse (set to 0 for SPGR)

if seqType<1, fail('seqType must be 1 (STFR) or 2 (SPGR)'), end
if seqType>2, fail('seqType must be 1 (STFR) or 2 (SPGR)'), end

[nx,ny,nz]=size(b0map);             % x, y, and z spatial dimensions
FOV=2*abs(simuRange.x(1));          % FOV (cm)
dt=4e-3;                            % dwell time (ms)
T1=Inf; T2=Inf;                     % set T1 and T2 to infinity to ignore relaxation effects
arg.targetf=25;                     % default local tracking bandwidth in Hz
arg.ktype=9;                        % default spectral-spatial excitation trajectory type
arg.undersamp=4;                    % default full spectral bandwidth
arg.spec_band=0;                    % default full spectral bandwidth
arg.Trf=3e-3;                       % default RF pulse length in sec
arg.sens=ones(1,nx,ny,nz);          % default sensitivity maps
arg.doplot=0;                       % default plotting
arg=vararg_pair(arg, varargin);     % reads in user-provided values
% if no excitation gradients, no local tracking bandwidth for
% spectral-spatial pulse
if arg.ktype==0
    arg.targetf=0;
end
% sets spectral bandwidth for purely spectral pulse
if arg.spec_band==1
    bwl=round(1e3/Tread)            % empirical spectral bandwidth limit
    fspec=-bwl/2:bwl/2;             % bandwidth window +/-bwl/2
else                                % use full spectral range of fieldmap
    fspec=0;
end

%% Design tip-down pulse

tstart1=tic;
theta_est=-2*pi*b0map*(Tread)/1e3;  % initial/nominal phase pattern
theta_est(~roi)=0;                  % windows out don't care region
magnom=ones(size(theta_est))*sin(nom_fa*pi/180);    % initial/nominal magnitude pattern
d_td_est=magnom.*exp(0.5*-1i*theta_est);            % initial/nominal design for prewinding pulse

if pulseType ==1                 	% spectral pulse design
   [btd,gxtd,gytd,gztd, mtd_approx, ftd, dtd, Atd, Wtd] = spectralRF(d_td_est,b0map,roi,arg.Trf,0,zslice,fspec);
   save Aspec.mat Atd Wtd;
elseif pulseType == 2           	% spectral-spatial pulse design
    [btd,gxtd,gytd,gztd,mtd_approx,ftd,dtd, circ,kxtd,kytd,kztd,Atd,Wtd]= spectralspatialRF(d_td_est,b0map,roi,arg.Trf,FOV,FOV,arg.undersamp,arg.undersamp,zslice,arg.targetf*1e-3,arg.ktype);
    [b12, gx2, gy2, gz2, mhat12, f12, d12, circ12, kx2, ky2, kz2, A12, W12]= spectralspatialRF(d_td_est, b0map, roi, Trf, FOV, FOV, undersamp, undersamp, zslice, targetf*1e-3, ktype);
    save Aspecspat.mat Atd Wtd;                                                                  
else
end
tipdownDesignTime=toc(tstart1);

%% Simulate after free precession and design tip-up pulse (if STFR)

if seqType==1
    % RF/gradient out to free precession
    nread = round(Tread/dt);            % number of time samples in readout (free precession time)
    b1d = [btd; zeros(nread,1)];
    gxd = [gxtd; zeros(nread,1)];
    gyd = [gytd; zeros(nread,1)];
    gzd = [gztd; zeros(nread,1)];
    % Bloch simulation of magnetization after free precession
    [m_beforeTipup,mz_beforeTipup]=parallel_blochCim(0,b1d,gxd,gyd,gzd,arg.sens(:,:,:,zslice), simuRange.x,...
        simuRange.y,simuRange.z(zslice),dt*1e-3,b0map(:,:,zslice),roi(:,:,zslice),T1*1e-3,T2*1e-3);
    d_tu_est=m_beforeTipup;                  % target pattern for tipup pulse
    d_tu_est=repmat(d_tu_est,[1 1 nz]);      % fills across 3D space (we only care about z-th slice)
    % design tip-up pulse to match magnetization pattern after free precession
    tstart2=tic;
    if pulseType == 1                   % spectral pulse
        [btu,gxtu,gytu,gztu, mtu_approx, ftu, dtu, Atu, Wtu] = spectralRF(d_tu_est, -1*b0map,roi, arg.Trf, 0, zslice,fspec);
        save Aspec.mat Atu Wtu -append;
    elseif pulseType == 2;              % spectral-spatial pulse
        [btu,gxtu,gytu,gztu,mtu_approx,ftu,dtu, circ,kxtu,kytu,kztu,Atu,Wtu]= spectralspatialRF(d_tu_est,-1*b0map,roi,arg.Trf,FOV,FOV,arg.undersamp,arg.undersamp,zslice,arg.targetf*1e-3,arg.ktype);
        save Aspecspat.mat Atu Wtu -append;
    end
    tipupDesignTime=toc(tstart2);
    btu=-flipud(btu);                  % Negate and time-reverse tip-up RF pulse and gradients
    gxtu=-flipud(gxtu);
    gytu=-flipud(gytu);
    gztu=-flipud(gztu);
elseif seqType==2
    btu=zeros(length(btd),1);           % No tip-up pulse for SPGR
    gxtu=zeros(length(gxtd),1);
    gytu=zeros(length(gytd),1);
    gztu=zeros(length(gztd),1);
    tipupDesignTime=0;
    m_beforeTipup=0;
end

%% Full pulse simulations: 1) At TE (STFR and SPGR) 2) After tip-up (STFR only)

% At TE
d=d_td_est(:,:,zslice);        % 2D target pattern
necho=round(Tread/2/dt);       % number of time samples at echo time (halfway through free precession)
b1e = [btd; zeros(necho,1)];   % simulate the excitation pattern at the echo time 
gxe = [-gxtd; zeros(necho,1)];        
gye = [-gytd; zeros(necho,1)];
gze = [-gztd; zeros(necho,1)];

m_atEcho=parallel_blochCim(0,b1e,gxe,gye,gze,arg.sens(:,:,:,zslice),simuRange.x,...
    simuRange.y,simuRange.z(zslice),dt*1e-3,b0map(:,:,zslice),roi(:,:,zslice),T1*1e-3,T2*1e-3);

if arg.doplot==1
    figure; 
    subplot(121); im(simuRange.x,simuRange.y,abs(m_atEcho)/sin(nom_fa*pi/180),[0 1]); 
    colormap(gca,'default'); xlabel('x (cm'); ylabel('y (cm)'); h=colorbar; 
    ylabel(h,'Normalized Magnitude'); title('Magnitude');
    subplot(122); im(simuRange.x,simuRange.y,angle(m_atEcho)*180/pi,[-180 180]); 
    colormap(gca,'hsv'); xlabel('x (cm'); ylabel('y (cm)'); h=colorbar; 
    ylabel(h,'Phase (\circ)'); title('Phase');
    suptitle('Magnetization at TE');
end
    
% After tip-up

if seqType==1
    b1f = [btd; zeros(nread,1); btu];% Full STFR RF pulse and gradient waveforms
    gxf = [gxtd; zeros(nread,1); gxtu];
    gyf = [gytd; zeros(nread,1); gytu];
    gzf = [gztd; zeros(nread,1); gztu];
    m_afterTipup=parallel_blochCim(0,b1f,gxf,gyf,gzf,arg.sens(:,:,:,zslice), simuRange.x,...
        simuRange.y,simuRange.z(zslice),dt*1e-3,b0map(:,:,zslice),roi(:,:,zslice),T1*1e-3,T2*1e-3);
    if arg.doplot==1
        figure; 
        subplot(121); im(simuRange.x,simuRange.y,abs(m_afterTipup)/sin(nom_fa*pi/180),[0 1]); 
        colormap(gca,'default'); xlabel('x (cm'); ylabel('y (cm)'); h=colorbar; 
        ylabel(h,'Normalized Magnitude'); title('Magnitude');
        subplot(122); im(simuRange.x,simuRange.y,angle(m_afterTipup)*180/pi,[-180 180]); 
        colormap(gca,'hsv'); xlabel('x (cm'); ylabel('y (cm)'); h=colorbar; 
        ylabel(h,'Phase (\circ)'); title('Phase');
        suptitle('Magnetization after Tip-up');
    end
else
    m_afterTipup=0;
end

return




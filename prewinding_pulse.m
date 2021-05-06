function [btd,gxtd,gytd,gztd,d,m_atEcho] = prewinding_pulse(pulseType, seqType, b0map, roi, simuRange, zslice, nom_fa, Tread, varargin)
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
% arg.Trf=6e-3;                       % (SLR) default RF pulse length in sec
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

nom_fa=90;                          % required flip angle in degrees
tstart1=tic;
theta_est=-2*pi*b0map*(Tread)/1e3;  % initial/nominal phase pattern
% theta_est=zeros(size(b0map)); % (SLR) assume no initial phase
theta_est(~roi)=0;                  % windows out don't care region
magnom=ones(size(theta_est))*sin(nom_fa*pi/180);    % initial/nominal magnitude pattern
d_td_est=magnom.*exp(0.5*-1i*theta_est);            % initial/nominal design for prewinding pulse

if pulseType ==1                 	% spectral pulse design
   [btd,gxtd,gytd,gztd, mtd_approx, ftd, dtd, Atd, Wtd] = spectralRF(d_td_est,b0map,roi,arg.Trf,0,zslice,fspec);
   save Aspec.mat Atd Wtd;
elseif pulseType == 2           	% spectral-spatial pulse design
    [btd,gxtd,gytd,gztd,mtd_approx,ftd,d, circ,kxtd,kytd,kztd,Atd,Wtd, time]= spectralspatialRF(d_td_est,b0map,roi,arg.Trf,FOV,FOV,arg.undersamp,arg.undersamp,zslice,arg.targetf*1e-3,arg.ktype);
    %[b12, gx2, gy2, gz2, mhat12, f12, d12, circ12, kx2, ky2, kz2, A12, W12]= spectralspatialRF(d_td_est, b0map, roi, Trf, FOV, FOV, undersamp, undersamp, zslice, targetf*1e-3, ktype);
    save Aspecspat.mat Atd Wtd;                                                                  
else
end
tipdownDesignTime=toc(tstart1);

%% Full pulse simulations: 1) At TE (STFR and SPGR) 2) After tip-up (STFR only)

% At TE
% d=d_td_est(:,:,zslice);        % 2D target pattern
% necho=round(Tread/2/dt);       % number of time samples at echo time (halfway through free precession)
necho=1;                       % (SLR) number of time samples at echo time (halfway through free precession)
b1e = [btd; zeros(necho,1)];   % simulate the excitation pattern at the echo time 
gxe = [gxtd; zeros(necho,1)];        
gye = [gytd; zeros(necho,1)];
gze = [gztd; zeros(necho,1)];

% m_atEcho=parallel_blochCim(0,b1e,gxe,gye,gze,arg.sens(:,:,:,zslice),simuRange.x,...
%    simuRange.y,simuRange.z(zslice),dt*1e-3,b0map(:,:,zslice),roi(:,:,zslice),T1*1e-3,T2*1e-3);
% =========================================================================
% SLR: simulate without the B0 map first

% simulate water magnetization
b0map_water_uniform = ones(size(b0map))*0;
m_atEcho=parallel_blochCim(0,b1e,gxe,gye,gze,arg.sens(:,:,:,zslice),simuRange.x,...
    simuRange.y,simuRange.z(zslice),dt*1e-3,b0map_water_uniform(:,:,zslice),roi(:,:,zslice),T1*1e-3,T2*1e-3);

% simulate fat by adding an offset to the B0 map
% b0map_fat = b0map - 3.3*127.74;
b0map_fat_uniform = ones(size(b0map))*(-3.3*127.74);
mFat_atEcho=parallel_blochCim(0,b1e,gxe,gye,gze,arg.sens(:,:,:,zslice),simuRange.x,...
    simuRange.y,simuRange.z(zslice),dt*1e-3,b0map_fat_uniform(:,:,zslice),roi(:,:,zslice),T1*1e-3,T2*1e-3);;

if arg.doplot==1
    figure; 
    subplot(121); im(simuRange.x,simuRange.y,abs(m_atEcho)/sind(nom_fa),[0 1]); 
    colormap(gca,'default'); xlabel('x (cm)'); ylabel('y (cm)'); h=colorbar; 
    ylabel(h,'Normalized Magnitude'); title('Magnitude');
    subplot(122); im(simuRange.x,simuRange.y,angle(m_atEcho)*180/pi,[-180 180]); 
    colormap(gca,'hsv'); xlabel('x (cm)'); ylabel('y (cm)'); h=colorbar; 
    ylabel(h,'Phase (\circ)'); title('Phase');
    sgtitle('Water magnetization after 1 time sample');
    
    figure; 
    subplot(121); im(simuRange.x,simuRange.y,abs(mFat_atEcho)/sind(nom_fa),[0 1]); 
    colormap(gca,'default'); xlabel('x (cm)'); ylabel('y (cm)'); h=colorbar; 
    ylabel(h,'Normalized Magnitude'); title('Magnitude');
    subplot(122); im(simuRange.x,simuRange.y,angle(mFat_atEcho)*180/pi,[-180 180]); 
    colormap(gca,'hsv'); xlabel('x (cm)'); ylabel('y (cm)'); h=colorbar; 
    ylabel(h,'Phase (\circ)'); title('Phase');
    sgtitle('Fat magnetization after 1 time sample');
    
    % Display the magnetization from small-tip angle approximation
    as(permute(abs(mtd_approx)/sind(nom_fa), [2 3 1]));
    as(permute(Wtd, [2 3 1]));
    as(permute(d, [2 3 1]));
    
    % RF magnitude, phase and gradients
    figure; 
    subplot(311); plot(time, abs(btd));
    xlabel('Time (ms)'); ylabel('B1 amplitude (Gauss)'); 
%     title('RF amplitude');
    subplot(312); plot(time, angle(btd)*180/pi);
    xlabel('Time (ms)'); ylabel('B1 phase (\circ)'); 
%     title('RF phase');
    subplot(313); plot(time, gxtd, 'r', time, gytd, 'g', time, gztd);
    xlabel('Time (ms)'); ylabel('Gradients (G/cm)');
    legend('G_x','G_y','G_z');
%     title('RF phase');
    sgtitle('RF pulse design');
    
    % k-space trajectory
    figure;
    plot3(kxtd, kytd, kztd); grid on
    xlabel('k_x (cm^{-1})'); ylabel('k_y (cm^{-1})'); zlabel('k_z (cm^{-1})');
    sgtitle('Pulse k-space trajectory');
end

% =========================================================================

return




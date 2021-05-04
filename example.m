%% Example of Prewinding RF Pulse Design
% runs example of spectral and spectral-spatial prewinding pulse design and simulations
%
% based on S. N. Williams, J-F. Nielsen, J.A. Fessler, and D.C. Noll,
% "Design of spectral-spatial phase prewinding pulses and their use in 
% small-tip fast recovery steady-state imaging", Magnetic Resonance in Medicine, 2017 DOI 10.1002/mrm.26794
%
% requires additional packages:
% Michigan IRT                          - Jeff Fessler's Image Recon Toolbox, http://web.eecs.umich.edu/~fessler/code/index.html
% Variable Density Spiral Design        - Brian Hargreaves' VDS Functions, http://www-mrsrl.stanford.edu/~brian/vdspiral/
% PTx Bloch Simulator                   - Hao Sun's Mex Bloch Simulator, http://www-personal.umich.edu/~sunhao/
%
% calls the following functions:
% prewinding_pulse.m                    - designs prewinding pulse (spectral or spectral-spatial)
%       -calls: spectralRF.m            - spectral pulse design
%                   qpwls_pcg1.m        - Fessler's penalized conjugate gradient function from Michigan IRT
%                   powerit.m           - Power iteration for estimating larget Eigenvalue
%                   ir_fista_syd.m      - FISTA algorithm for pulse design with peak/power constraints, adapted from ir_fista.m
%               spectralspatialRF.m     - spectral-spatial pulse design
%                   excite_vds_grad.m   - designs various variable density spiral trajectories for excitation
%                   g2k.m               - converts from gradient trajectory to k-space trajectory
%                   qpwls_pcg1.m        - Fessler's penalized conjugate gradient function from Michigan IRT
%                   powerit.m           - Power iteration for estimating larget Eigenvalue
%                   ir_fista_syd.m      - FISTA algorithm for pulse design with peak/power constraints, adapted from ir_fista.m
%               parallel_blochCim.m     - runs Bloch simulation with PTx capacity, written by H.Sun 2013 adapted from B. Hargreaves, blochCim.c (mex) 
% bloch_4D.m                            - runs Bloch simulation of magnetization over x,y,z as well as multiple intravoxel frequencies
%       -calls: weight_matrix.m         - designs weighting matrix based on target intravoxel range and shape
%                   gauss_weight.m      - creates vector of Gaussian weights for particular sigma standard deviation
%               parallel_Bloch4D.m      - 4D implementation of parallel_blochCim.m for multiple intravoxel frequencies
% combine_3D.m                          - combines multi-frequency Bloch simulation for single 2D slice
%       -calls: freq_vector.m           - builds frequency vector for particular range and sampling rate
%               gauss_weight.m          - creates vector of Gaussian weights for particular sigma standard deviation
% perform_metrics.m                     - computes excitation NRMSE, phase RMSE, magnitude NRMSE, mean and standard deviation of magnitude for particular pulse and target pattern
% bloch_1spat1freq.m                    - Bloch simulation for examining full-range of intravoxel frequencies and one spatial dimension
%       -calls: weight_matrix.m         - designs weighting matrix based on target intravoxel range and shape
%                   gauss_weight.m      - creates vector of Gaussian weights for particular sigma standard deviation
%               parallel_Bloch4D.m      - 4D implementation of parallel_blochCim.m for multiple intravoxel frequencies
% plot_human_data.m                     - plots experimental and simulation human data
% plot_spatVfreq.m                      - plots frequency simulation vs 1 spatial dimension plots 
%
load PrewindingPulse.mat                                % loads in pre-computed field map and some acquisition/design details
pulseType=1;                                            % RF pulse type for design, 1: purely spectral 2: spectral-spatial
seqType=2;                                              % pulse sequence, 1: STFR 2: SPGR 
%% Pulse Design
% spectral
[btd_spec,gxtd_spec,gytd_spec,gztd_spec,...         % designs tip-down (and tip-up) prewinding pulse set
    btu_spec,gxtu_spec,gytu_spec,gztu_spec,...
    d,m_atEcho_spec,m_beforeTipup_spec,m_afterTipup_spec]=prewinding_pulse(pulseType,...
    seqType, b0map, roi, simuRange,zslice,nom_fa, Tread, 'doplot',1);
% spectral-spatial
% designs tip-down (and tip-up) prewinding pulse set
[btd_specspat,gxtd_specspat,gytd_specspat,gztd_specspat, btu_specspat,gxtu_specspat,gytu_specspat,gztu_specspat, d,m_atEcho_specspat,m_beforeTipup_specspat,m_afterTipup_specspat]=prewinding_pulse((pulseType+1), seqType, b0map, roi, simuRange,zslice,nom_fa, Tread, 'doplot',1);
% here, we want uniform magnetization and zero phase at TE
% we also want zero magnetization after tip-up (don't care about phase)
%% Bloch Simulation/Performance Metrics
fsamp=10;                                               % frequency sampling rate for intravoxel simulation, typically 10 Hz
% 4D Bloch simulation across (x,y,z, and frequency space) --3D here for 2D design
% spectral
[mBloch_spec,fBloch_spec,b0map_tp_spec]=bloch_4D(btd_spec,...
    gxtd_spec,gytd_spec,gztd_spec,b0map,roi,simuRange,Tread,...
    'zslice',zslice,'fsamp',fsamp);
% spectral-spatial
[mBloch_specspat,fBloch_specspat,b0map_tp_specspat]=bloch_4D(btd_specspat,...
   gxtd_specspat,gytd_specspat,gztd_specspat,b0map,roi,simuRange,Tread,...
   'zslice',zslice,'fsamp',fsamp);
% hard pulse
[mBloch_hard,fBloch_hard,b0map_tp_hard]=bloch_4D(b1_hard,...
    gxtd_hard,gytd_hard,gztd_hard,b0map,roi,simuRange,Tread,...
    'zslice',zslice,'fsamp',fsamp);

% combines multi-dimensional Bloch simulation across intravoxel frequency space
% spectral
[meanFBloch_spec,abs_meanFBloch_spec,angle_meanFBloch_spec]=...
    combine_3D(mBloch_spec,fBloch_spec,b0map_tp_spec,fsamp);
% spectral-spatial
[meanFBloch_specspat,abs_meanFBloch_specspat,angle_meanFBloch_specspat]=...
    combine_3D(mBloch_specspat,fBloch_specspat,b0map_tp_specspat,fsamp);
% hard pulse
[meanFBloch_hard,abs_meanFBloch_hard,angle_meanFBloch_hard]=...
    combine_3D(mBloch_hard,fBloch_hard,b0map_tp_hard,fsamp);

% computes performance metrics
% spectral
[mean_mxy_spec,perstd_mxy_spec,sNRMSE_mxy_spec,magNRMSE_mxy_spec,...
    phsRMSE_mxy_spec]=perform_metrics(meanFBloch_spec,...
    abs_meanFBloch_spec,angle_meanFBloch_spec,d,nom_fa,...
    roi(:,:,zslice),'doplot',1);
% spectral-spatial
[mean_mxy_specspat,perstd_mxy_specspat,sNRMSE_mxy_specspat,magNRMSE_mxy_specspat,...
    phsRMSE_mxy_specspat]=perform_metrics(meanFBloch_specspat,...
    abs_meanFBloch_specspat,angle_meanFBloch_specspat,d,nom_fa,...
    roi(:,:,zslice),'doplot',1);
% hard pulse
[mean_mxy_hard,perstd_mxy_hard,sNRMSE_mxy_hard,magNRMSE_mxy_hard,...
    phsRMSE_mxy_hard]=perform_metrics(meanFBloch_hard,...
    abs_meanFBloch_hard,angle_meanFBloch_hard,d,nom_fa,...
    roi(:,:,zslice),'doplot',1);

% spatial vs all frequency simulations (from Supporting Info)
fsamp_fine=5;
% hard pulse
[mBlochF_hard,f_hard,b0map_center]=bloch_1spat1freq(b1_hard,gxtd_hard,...
    gytd_hard,gztd_hard,1,b0map,roi,simuRange,Tread,'zslice',zslice,'fsamp',fsamp_fine);
% spectral pulse
[mBlochF_spec,f_spec,b0map_center]=bloch_1spat1freq(btd_spec,gxtd_spec,...
    gytd_spec,gztd_spec,0,b0map,roi,simuRange,Tread,'zslice',zslice,'fsamp',fsamp_fine);
% spectral-spatial pulse
[mBlochF_specspat, f_specspat]=bloch_1spat1freq(btd_specspat,gxtd_specspat,...
    gytd_specspat,gztd_specspat,0,b0map,roi,simuRange,Tread,'zslice',zslice,'fsamp',fsamp_fine);
%% Plots

% simulation/experimental plots
[sNRMSE_space_hard, sNRMSE_space_spec, sNRMSE_space_specspat, ...
    phsRMSE_space_hard, phsRMSE_space_spec, phsRMSE_space_specspat]...
    =plot_human_data(meanFBloch_spec, meanFBloch_specspat, ...
    meanFBloch_hard, roi(:,:,zslice), d, simuRange.x, simuRange.y, nom_fa, im_spec,im_specspat, ...
    relphase_spec,relphase_specspat);

x1=62; y1=75;                        % coordinates for line profiles
xmi=30; xma=88; ymi=23; yma=113;     % range for plotting
targetf=mean(b0map_tp_spec(roi(:,:,zslice)));% target intravoxel bandwidth in Hz                          

% plots one spatial dimension vs frequency
% hard pulse
[magE_hard,phaseE_hard,complexE_hard,fsim_hard]=plot_spatVfreq(...
    simuRange.x,simuRange.y,xmi,xma,ymi,yma,x1,y1,fliplr(b0map_center),...
    fsamp_fine,nom_fa,mBlochF_hard,abs_meanFBloch_hard, ...
    angle_meanFBloch_hard, f_hard, targetf);
% spectral pulse
[magE_spec,phaseE_spec,complexE_spec,fsim_spec]=plot_spatVfreq(...
    simuRange.x,simuRange.y,xmi,xma,ymi,yma,x1,y1,fliplr(b0map_center),...
    fsamp_fine,nom_fa,mBlochF_spec,abs_meanFBloch_spec, ...
    angle_meanFBloch_spec, f_spec, targetf);
% spectral-spatial pulse
[magE_specspat,phaseE_specspat,complexE_specspat,fsim_specspat]=plot_spatVfreq(...
    simuRange.x,simuRange.y,xmi,xma,ymi,yma,x1,y1,fliplr(b0map_center),...
    fsamp_fine,nom_fa,mBlochF_specspat,abs_meanFBloch_specspat, ...
    angle_meanFBloch_specspat, f_specspat, targetf);

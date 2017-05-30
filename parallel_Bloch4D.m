function [Mxy_full, Mz_full,f]=parallel_Bloch4D(Minit,pulses,gx,gy,gz,sens,x,y,z,dt,b0map,roi, T1, T2,minf,maxf,samp_f,mode)
% function [Mxy_full, Mz_full,f]=parallel_Bloch4D(Minit,pulses,gx,gy,gz,sens,x,y,z,dt,b0map,roi, T1, T2,minf,maxf,samp_f,mode)
% parallel Bloch simulation for 4D (3 spatial, one spectral)
% loops through Hao Sun's Parallel Bloch Simulator (2012) adapted from
% Brian Hargreave's website
% S. Williams, University of Michigan, 03.02.2016

% Inputs:
% Minit         [1]             0 --> all magnetization in Mz [0, 0, 1] (unless known different!)
% pulses        [npts nc]       complex valued RF pulse
% gx            [npts 1]        excitation gradients in X
% gy            [npts 1]        excitation gradients in Y
% gz            [npts 1]        excitation gradients in Z
% sens          [nc nx ny nz]   b1 maps
% x             [nx]            spatial samples over FOV in X
% y             [ny]            spatial samples over FOV in Y
% z             [nz]            spatial samples over FOV in Z
% dt            [1]             dwell time sampling RF pulse (ms)
% b0map         [nx ny nz]      off-resonance fieldmap (Hz)
% roi           [nx ny nz]      excitation map (logical)
% T1            [1]             T1 relaxation time (ms)
% T2            [1]             T2 decay time (ms)
% minf          [1]             smallest frequency to begin simulation at (Hz)
% maxf          [1]             largest frequency to end simulation at (Hz)
% samp_f        [1]             spectral sampling frequency (Hz)
% mode          [1]             0: just 1 TR 1: steady-state (careful, expensive!)

% Outputs:
% Mxy_full      [nx ny nz nf]  simulated transverse magnetization across all frequencies in spectral design
% Mz_full       [nx ny nz nf]  simulated longitudinal magnetization across all frequencies in spectral design
% f             [nf]           vector of frequency shift range used in simulation (Hz)

    f=freq_vector(minf,maxf,samp_f);
    nf=length(f);
    Mxy_full=zeros([size(b0map) nf]);       % placeholder for transver/longitudinal magnetizations
    Mz_full=zeros([size(b0map) nf]);
    for i=1:nf
        b0tmp=b0map+f(i);                   % shifted fieldmap
        [mxy,mz] = parallel_blochCim(Minit,pulses,gx,gy,gz,sens,x,y,z,...
            dt*1e-3,b0tmp,roi,T1*1e-3,T2*1e-3,mode);
        if ndims(mxy)==3
            Mxy_full(:,:,:,i)=mxy;              % transverse magnetization at f(i)
            Mz_full(:,:,:,i)=mz;                % longitudinal magnetization at f(i)
        elseif ndims(mxy)<3
            Mxy_full(:,:,i)=mxy;              % transverse magnetization at f(i)
            Mz_full(:,:,i)=mz;                % longitudinal magnetization at f(i)
        end
        sprintf('niter: %2.1f out of %2.1f',i,nf+1)
    end

end
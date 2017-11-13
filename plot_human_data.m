function [sNRMSE_space_hard, sNRMSE_space_spec, sNRMSE_space_specspat, phsRMSE_space_hard, phsRMSE_space_spec, phsRMSE_space_specspat] = plot_human_data(meanFBloch_spec, meanFBloch_specspat, meanFBloch_hard, roi, d, x, y, nom_fa, im_spec,im_specspat, relphase_spec,relphase_specspat)
% function [sNRMSE_space_hard, sNRMSE_space_spec, sNRMSE_space_specspat, phsRMSE_space_hard, phsRMSE_space_spec, phsRMSE_space_specspat] nplot_human_data(meanFBloch_spec, meanFBloch_specspat, meanFBloch_hard, roi, d, x, y, nom_fa, im_spec,im_specspat, relphase_spec,relphase_specspat)
% reproduces simulation and experimental images for human data from
% Sydney Williams, University of Michigan 2017
%
%
% INPUTS:
% (required)
% meanFBloch_spec       [nx ny]   	- simulated multi-frequency Bloch simulation for spectral pulse
% meanFBloch_specspat   [nx ny]  	- simulated multi-frequency Bloch simulation for spectral-spatial pulse
% meanFBloch_hard       [nx ny]    	- simulated multi-frequency Bloch simulation for hard pulse
% roi                   [nx ny]    	- logical ROI of 2D slice
% d                     [nx ny]    	- target pattern for 2D slice
% x                     [nx]       	- x vector in cm
% y                     [ny]       	- y vector in cm
% nom_fa                [1]        	- target flip angle in degrees
% im_spec               [nxh nyh]   - high resolution experimental purely spectral image
% im_specspat           [nxh nyh]   - high resolution experimental spectral-spatial image
% relphase_spec         [nxh nyh]   - high resolution experimental purely spectral relative phase image (common phase removed)
% relphase_specspat     [nxh nyh]   - high resolution experimental spectral-spatial relative phase image (common phase removed)
%
% OUTPUTS:
% sNRMSE_space_hard     [nx ny]    	- excitation NRMSE as a function of space for hard pulse
% sNRMSE_space_spec    	[nx ny]   	- excitation NRMSE as a function of space for spectral pulse
% sNRMSE_space_specspat [nx ny]  	- excitation NRMSE as a function of space for spectral-spatial pulse
% phsRMSE_space_hard  	[nx ny]    	- phase RMSE as a function of space for hard pulse
% phsRMSE_space_spec 	[nx ny]   	- phase RMSE as a function of space for spectral pulse
% phsRMSE_space_specspat[nx ny]  	- phase RMSE as a function of space for spectral-spatial pulse
%

set(0,'defaultAxesFontName', 'Times')
set(0,'defaultTextFontSize', 13)
set(0,'defaultTextFontName', 'Times')
set(0,'defaultAxesFontSize', 13)

% plots simulation figures
figure; 
subplot(221); im(x,y,(meanFBloch_spec)/sin(nom_fa*pi/180),[0 1]); colormap(gca,'gray');
    ylabel('y (cm)'); title('(a)'); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'}); 
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'});  
    set(gca,'Position', [0.3130 0.5621 0.1825 0.3625]);
    ylabh = get(gca,'YLabel'); set(ylabh,'Position',[-18 -20 1]);
subplot(222); im(x,y,(meanFBloch_specspat)/sin(nom_fa*pi/180),'cbar',[0 1]); colormap(gca,'gray');
    title('(b)'); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'});
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); 
    c=colorbar; ylabel(c,{'Normalized', 'magnitude'}); set(c,'ytick',[0 1],'yTickLabel',{'0','1'})
    set(gca,'Position', [0.5513 0.5621 0.1825 0.3625]);
subplot(223); im(x,y,(angle(meanFBloch_spec)*180/pi),[-180 180]); colormap(gca,'hsv');
    xlabel('x (cm)'); title('');set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'}); xlabel('x (cm)');
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'});
    set(gca,'Position', [0.3130 0.0980 0.1825 0.3625]);
    xlabh=get(gca,'XLabel'); set(xlabh,'Position',[16 -16 1]);
subplot(224); im(x,y,(angle(meanFBloch_specspat)*180/pi),'cbar',[-180 180]); colormap(gca,'hsv');
    title(''); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'}); 
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); 
    c=colorbar; ylabel(c,'Phase (\circ)'); set(c,'ytick',[-180 0 180],'yTickLabel',{'180','0','180'});
    set(gca,'Position', [0.5513 0.0980 0.1825 0.3625]);
suptitle('Bloch Simulated Images');
% plots experimental images
[nxh,nyh]=size(im_spec);                    % computes x-y vectors for high(er) resolution experimental images
FOV=24; 
xh=[-nxh/2:(nxh/2-1)]*(FOV/nxh);
yh=[-nyh/2:(nyh/2-1)]*(FOV/nyh);

figure;
subplot(221); im(xh,yh,fliplr(im_spec),[0 0.09]); colormap(gca,'gray');
    ylabel('y (cm)'); title('(a)'); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'}); 
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); 
    c=colorbar; ylabel(c,'Magnitude (a.u.)'); set(c,'ytick',[0 0.15],'yTickLabel',{'0','0.09'})
    set(gca,'Position', [0.2350 0.5621 0.1825 0.3625]);
    ylabh = get(gca,'YLabel'); set(ylabh,'Position',[-18 -20 1]);
subplot(222); im(xh,yh,fliplr(im_specspat),'cbar',[0 0.15]); colormap(gca,'gray');
    title('(b)'); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'});
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); 
    c=colorbar; ylabel(c,'Magnitude (a.u.)');  set(c,'ytick',[0 0.15],'yTickLabel',{'0','0.15'});
    set(gca,'Position', [0.5900 0.5621 0.1825 0.3625]);
subplot(223); im(xh,yh,fliplr(relphase_spec),[-180 180]); colormap(gca,'hsv');
    xlabel('x (cm)'); title('');set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'}); xlabel('x (cm)');
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); 
    c=colorbar; ylabel(c,'Phase (\circ)'); set(c,'ytick',[-180 0 180],'yTickLabel',{'180','0','180'});
    set(gca,'Position', [0.2350 0.0980 0.1825 0.3625]);
    xlabh=get(gca,'XLabel'); set(xlabh,'Position',[25 -16 1]);
subplot(224); im(xh,yh,fliplr(relphase_specspat),'cbar',[-180 180]); colormap(gca,'hsv');
    title(''); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'}); 
    set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); 
    c=colorbar; ylabel(c,'Phase (\circ)'); set(c,'ytick',[-180 0 180],'yTickLabel',{'180','0','180'});
    set(gca,'Position', [0.5900 0.0980 0.1825 0.3625]);
suptitle('Experimental Images');  

% excitation NRMSE and Phase RMSE over spatial locations
angle_meanFBloch_hard=angle(meanFBloch_hard)*180/pi;
angle_meanFBloch_spec=angle(meanFBloch_spec)*180/pi;
angle_meanFBloch_specspat=angle(meanFBloch_specspat)*180/pi;

sNRMSE_space_hard=sqrt((abs(meanFBloch_hard.*roi-abs(d.*roi)).^2))...
    /(sin(nom_fa*pi/180));                                      % excitation NRMSE as function of space  
sNRMSE_space_spec=sqrt((abs(meanFBloch_spec.*roi-abs(d.*roi)).^2))...
    /(sin(nom_fa*pi/180));                                      
sNRMSE_space_specspat=sqrt((abs((meanFBloch_specspat).*roi-abs(d.*roi)).^2))...
    /(sin(nom_fa*pi/180));                                      
sNRMSE_hard=sqrt(mean(abs(meanFBloch_hard(roi)-abs(d(roi))).^2))...
        /(sin(nom_fa*pi/180));                                  % excitation NRMSE (single metric)
sNRMSE_spec=sqrt(mean(abs(meanFBloch_spec(roi)-abs(d(roi))).^2))...
        /(sin(nom_fa*pi/180));                              
sNRMSE_specspat=sqrt(mean(abs(meanFBloch_specspat(roi)-abs(d(roi))).^2))...
        /(sin(nom_fa*pi/180));                                

phsRMSE_space_hard=sqrt((angle_meanFBloch_hard.*roi).^2);       % phase RMSE as function of space       
phsRMSE_space_spec=sqrt((angle_meanFBloch_spec.*roi).^2); 
phsRMSE_space_specspat=sqrt(((angle_meanFBloch_specspat).*roi).^2);
phsRMSE_hard=sqrt(mean(angle_meanFBloch_hard(roi).^2));        	% phase RMSE (single metric)      
phsRMSE_spec=sqrt(mean(angle_meanFBloch_spec(roi).^2));
phsRMSE_specspat=sqrt(mean(angle_meanFBloch_specspat(roi).^2));

figure; % all three pulses
% hard pulse
subplot(231); im(x,y,(sNRMSE_space_hard),[0 1]); colormap(gca,'default'); 
set(gca,'Position',[0.1000	0.5450 0.2150	0.3350])
xticks([-12 0 11.9]); xticklabels({'-12','0','12',}); 
ylabel('y (cm)'); yticks([-12 0 11.9]); yticklabels({'-12','0','12',}); 
ylabh = get(gca,'YLabel'); set(ylabh,'Position',[-17 -16 1]);
title(sprintf('Hard Pulse \n E.NRMSE=%3.2f',sNRMSE_hard));
subplot(234); im(x,y,(phsRMSE_space_hard),[0 90]); colormap(gca,'default'); 
set(gca,'Position',[0.1000	0.0950 0.2150	0.3350])
xlabel('x (cm)'); xticks([-12 0 11.9]); xticklabels({'-12','0','12',}); 
yticks([-12 0 11.9]); yticklabels({'-12','0','12',}); 
xlabh = get(gca,'XLabel'); set(xlabh,'Position',[33 -15 1]);
title(sprintf('Hard Pulse \n Phase RMSE=%3.1f',phsRMSE_hard));
% purely spectral pulse
subplot(232); im(x,y,(sNRMSE_space_spec),[0 1]); colormap(gca,'default');
set(gca,'Position',[0.3750	0.5450 0.2150	0.3350])
xticks([-12 0 11.9]); xticklabels({'-12','0','12'}); 
yticks([-12 0 11.9]); yticklabels({'-12','0','12'});
title(sprintf('Purely Spectral Pulse \n E. NRMSE=%3.2f',sNRMSE_spec));
subplot(235); im(x,y,((phsRMSE_space_spec)),[0 90]); colormap(gca,'default');
set(gca,'Position',[0.3750	0.0950 0.2150	0.3350])
xticks([-12 0 11.9]); xticklabels({'-12','0','12'}); 
yticks([-12 0 11.9]); yticklabels({'-12','0','12'}); 
title(sprintf('Purely Spectral Pulse \n Phase RMSE=%3.1f',phsRMSE_spec));
% spectral-spatial pulse
subplot(233); im(x,y,((sNRMSE_space_specspat)),[0 1]); h=colorbar; colormap(gca,'default');
set(gca,'Position',[0.6500	0.5450 0.2150	0.3350])
xticks([-12 0 11.9]); xticklabels({'-12','0','12'}); 
yticks([-12 0 11.9]); yticklabels({'-12','0','12'}); 
set(h,'YTick',[0 0.5 1]);
title(sprintf('Spectral-Spatial Pulse \n E. NRMSE=%3.2f',sNRMSE_specspat));
subplot(236); im(x,y,((phsRMSE_space_specspat)),[0 90]); h=colorbar; colormap(gca,'default');
set(gca,'Position',[0.6500	0.0950 0.2150	0.3350])
xticks([-12 0 11.9]); xticklabels({'-12','0','12'}); 
yticks([-12 0 11.9]); yticklabels({'-12','0','12'}); 
set(h,'YTick',[0 45 90]);
title(sprintf('Spectral-Spatial Pulse \n Phase RMSE=%3.1f',phsRMSE_specspat));

return
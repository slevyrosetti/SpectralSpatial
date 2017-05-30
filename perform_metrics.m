function [mean_mxy, perstd_mxy, sNRMSE_mxy, magNRMSE_mxy, phsRMSE_mxy] = perform_metrics(meanFBloch, abs_meanFBloch, angle_meanFBloch, d, nom_fa, roiz, varargin)
% function [mean_mxy, perstd_mxy, sNRMSE_mxy, magNRMSE_mxy, phsRMSE_mxy] = perform_metrics(meanFBloch, abs_meanFBloch, angle_meanFBloch, d, nom_fa, roiz, varargin)
% Takes in Bloch simulated magnetization and target pattern, computes
% various performance metrics and has option for plot
%  Sydney Williams, University of Michigan 2016
%
% INPUTS:
% (required)
% meanFBloch       [nx ny]       - mean complex transverse magnetization from Bloch
%                                  simulation, potentially 4D simulation combined 
%                                  with combine_4D
% abs_meanFBloch   [nx ny]       - mean transverse magnetization magnitude from Bloch
%                                  simulation, potentially 4D simulation combined 
%                                  with combine_4D
% angle_meanFBloch [nx ny]       - mean transverse magnetization phase from Bloch
%                                  simulation, potentially 4D simulation combined 
%                                  with combine_4D
% d                [nx ny]       - target 2D magnetization pattern
% nom_fa           [1]           - target flip angle in degrees
% roiz             [nx ny]       - 2D ROI for slice
%
% (optional)
% 'FOV'            [1]           - FOV in cm
% 'doplot'         [0|1]         - set to one to add plot 
%
% OUTPUTS:
% mean_mxy         [1]           - mean magnetization magnitude
% perstd_mxy       [1]           - percent standard deviation of magnitude
% sNRMSE_mxy       [1]           - excitation NRMSE as defined in specspat paper
% magNRMSE_mxy     [1]           - magnitude NRMSE as defined in specspat paper
% phsRMSE_mxy      [1]           - phase RMSE as defined in specspat paper

    % setup parameters
    arg.FOV=24;                             % default FOV in cm (assumes same in x,y dimensions)
    arg.doplot=0;                           % default is to not plot
    arg = vararg_pair(arg, varargin);
    [nx,ny]=size(meanFBloch);               % size of 2D slice
    x=[-nx/2:(nx/2-1)]*(arg.FOV/nx);        % x sample locations (cm)
    y=[-ny/2:(ny/2-1)]*(arg.FOV/ny);        % y sample locations (cm)

    % computes performance metrics
    mean_mxy=mean(abs_meanFBloch(roiz));    % mean magnetization magnitude
    perstd_mxy=std(abs_meanFBloch(roiz))... 
        /(mean_mxy)*100;                    % percent standard devation of magnitude
    sNRMSE_mxy=sqrt(mean(abs(meanFBloch(roiz)-abs(d(roiz))).^2))...
        /(sin(nom_fa*pi/180));              % excitation NRMSE
    magNRMSE_mxy=sqrt(mean((abs_meanFBloch(roiz)-abs(d(roiz))).^2))...
        /sin(nom_fa*pi/180);                % magnitude NRMSE      
    phsRMSE_mxy=sqrt(mean(angle_meanFBloch(roiz).^2));% phase RMSE                           
    
    % optional plots
    if arg.doplot
        set(0,'defaultTextFontSize', 12);
        set(0,'defaultAxesFontSize', 12);
        set(0,'defaultAxesFontName', 'Times');
       % set(0,'defaultFontName', 'Times');
        figure; 
        subplot(121); im(x,y,fliplr(abs_meanFBloch)/sin(nom_fa*pi/180),'cbar',[0 1],...
            sprintf('Excitation \n NRMSE=%1.2f',sNRMSE_mxy)); colormap(gca,'gray'); 
            ylabel('y (cm)'); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'}); 
            set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); xlabel('x (cm)');
            set(gca,'Position',[0.1050    0.1010    0.3050    0.7500]);
            h=colorbar; ylabel(h,'Normalized Magnitude','FontSize',12,'FontName','Times'); 
            set(h,'ytick',[0 0.5 1],'yTickLabel',{'0','0.5','1.0'},'FontSize',12);
        subplot(122); im(x,y,fliplr(angle_meanFBloch),[-180 180],'cbar',...
            sprintf('Bloch \n phase RMSE=%2.1f%c',phsRMSE_mxy,char(176))); colormap(gca,'hsv');
            ylabel('y (cm)'); set(gca,'xtick',[-12 0 11.8],'xTickLabel',{'-12', '0','12'});
            set(gca,'ytick',[-12 0 11.8],'yTickLabel',{'-12', '0','12'}); xlabel('x (cm)');
            set(gca,'Position',[0.5700    0.1015    0.3050    0.7500]);
            h=colorbar; ylabel(h,'Phase (\circ)','FontSize',12,'FontName','Times'); 
            set(h,'ytick',[-180 0 180],'yTickLabel',{'-180','0','180'},'FontSize',12);
        suptitle(sprintf('Mean Magnitude=%1.2f, Percent St.Dev=%2.1f%c, \n Magnitude NRMSE=%1.2f',...
            mean_mxy, perstd_mxy, char(37), magNRMSE_mxy));
    end
end
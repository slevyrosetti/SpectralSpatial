function [magE, phaseE, complexE, fsim] = plot_spatVfreq(x,y,xmi,xma,ymi,yma,x1,y1,b0map_val,samp_f,nom_fa,mF,abs_meanFBloch, angle_meanFBloch, f, targetf)
% function plot_spatVfreq(x,y,xmi,xma,ymi,yma,x1,y1,b0map_val,samp_f,nom_fa,mF,abs_meanFBloch, angle_meanFBloch, f, targetf)
% plots "cool" frequency vs 1 spatial dimension plots from Williams et al.,MRM 2017
% Sydney Williams, University of Michigan 2017
%
% Inputs:
% x             [Nx]        -vector of x-locations (cm)
% y             [Ny]        -vector of x-locations (cm)
% xmi           [1]         -minimum plotting location index for x
% xma           [1]         -maximum plotting location index for x
% ymi           [1]         -minimum plotting location index for y
% yma           [1]         -maximum plotting location index for y
% x1            [1]         -location value of field map to span 1D profile across y
% x1            [1]         -location value of field map to span 1D profile across x
% b0map_val     [Nx Ny]     -masked field map of slice shifted to be centered at median frequency (Hz)
% samp_f        [1]         -frequency sampling rate (Hz)
% nom_fa        [1]         -design target flip angle (degrees)
% mF            [Nx Ny Nf]  -multiple frequency simulations
% f             [Nf]        -frequency vector used in mF
% abs_meanFBloch [Nx Ny]    -combined multifrequency Bloch simulation magnitude
% angle_meanFBloch [Nx Ny]  -combined multifrequency Bloch simulation magnitude
% targetf       [Nf]        -spectral-spatial design frequency bandwidth
%
% Outputs:
% magE          [Nx Ny Nf]  -relative magnitude error for each frequency ||m_xy|/sin(alpha)-1|
% phaseE        [Nx Ny Nf]  -relative phase error for each frequency |<m_xy|
% complexE      [Nx Ny Nf]  -relative complex error for each frequency |m_xy/sin(alpha)-1|
% fsim          [Nf]        -vector of frequencys f used


set(0,'defaultAxesFontName', 'Times')
set(0,'defaultAxesFontSize', 13)
set(0,'defaultTextFontName', 'Times')
set(0,'defaultAxesFontSize', 16)

magE=abs(abs(mF)/sin(nom_fa*pi/180)-1);
phaseE=abs(angle(mF)*180/pi);
complexE=abs(mF/sin(nom_fa*pi/180)-1);
fsim=freq_vector(-3*targetf,3*targetf,samp_f);
b0map_x=(5*round(b0map_val(xmi:xma,y1)/5)+targetf);
b0map_y=(5*round(b0map_val(x1,ymi:yma)/5)+targetf);
[nx,ny]=size(b0map_val);

figure; %orient tall; % relative error: magnitude, phase, and complex
% plots field map of slice
subplot(3,3,1);im(x(ymi:yma),y(ymi:yma),fliplr(b0map_val(ymi:yma,ymi:yma))); hold on;
%set(gca,'Position',[0.1300    0.7522    0.1578    0.1519]);
set(gca,'Position',[0.1300 0.7360 0.2000 0.2000]);
xticks([-7.5 0 10]); xticklabels({'-7.5','0','10'});
ylabel('y (cm)'); yticks([-7.5 0 10]); yticklabels({'-7.5','0','10'}); h=colorbar; 
colormap(gca,'jet'); ylabel(h,'Hz','fontsize',16); line([x(x1) x(x1)],[y(ymi) y(yma)],'LineWidth',2.0,'Color',[0 0 0]);
line([x(ymi) x(yma)],[y(y1) y(y1)],'LineWidth',2.0,'Color',[0 0 0]); hold off; title('Fieldmap With Line Profiles'); 
% plots simulated magnetization magnitude of slice
subplot(3,3,2); im(x(ymi:yma),y(ymi:yma),fliplr(abs_meanFBloch(ymi:yma,ymi:yma))/sin(nom_fa*pi/180),[0 1]); 
%set(gca,'Position',[0.4108    0.7522    0.1578    0.1519]);
set(gca,'Position',[0.4330 0.7360 0.2000 0.2000]);
colormap(gca,'default'); xlabel('x (cm)'); xticks([-7.5 0 10]); xticklabels({'-7.5','0','10'});
yticks([-7.5 0 10]); yticklabels({'-7.5','0','10'}); h=colorbar; 
ylabel(h,'Normalized Magnitude','fontsize',16); title('Simulated Magnitude'); 
% plots simulated magnetization phase of slice    
subplot(3,3,3); im(x(ymi:yma),y(ymi:yma),fliplr(angle_meanFBloch(ymi:yma,ymi:yma)),[-180 180]); 
%set(gca,'Position',[0.6916    0.7522    0.1578    0.1519]);
set(gca,'Position',[0.7360 0.7360 0.2000 0.2000]);
colormap(gca,'hsv'); xticks([-7.5 0 10]); xticklabels({'-7.5','0','10'});
yticks([-7.5 0 10]); yticklabels({'-10','0','7.5'}); h=colorbar; 
ylabel(h,'Phase (\circ)','fontsize',16); title('Simulated Phase'); 

% one y-location, complex error
subplot(334); im(x(xmi:xma),f,fliplr(squeeze(complexE((nx-xma):(nx-xmi),y1,:))),[0 3]); hold on;
%set(gca,'Position',[0.1308    0.4583    0.1498    0.1425]);
set(gca,'Position',[0.1300 0.3900    0.1250    0.2150]);
xticks([-5 0 5]); xticklabels({'-5','0','5'});
yticks([-115 0 100 205]); yticklabels({'-115','0','100','205'});
colormap(gca,'default'); ylabel('f (Hz)'); h=colorbar; ylabel(h,'$|m_{\mathrm{xy}}/\mathrm{sin}(\alpha)-1|$', 'interpreter', 'latex', 'fontsize', 16);
plot(x(xmi:xma),b0map_x,'LineWidth',2.0,'Color',[0 0 0]); plot(x(xmi:xma),b0map_x+targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);
plot(x(xmi:xma),b0map_x-targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);plot(x(xmi:xma),b0map_x+max(fsim),'r--','LineWidth',1.0);
plot(x(xmi:xma),b0map_x+min(fsim),'r--','LineWidth',1.0); hold off; title(sprintf('Rel. Complex Error at Y=%2.1f cm',y(y1))); 
%title(sprintf('$$Rel. Complex Error at Y=%2.1f cm$$',y(y1)),'interpreter','latex');  
% one y-location, magnitude error
subplot(335); im(x(xmi:xma),f,fliplr(squeeze(magE((nx-xma):(nx-xmi),y1,:))),[0 1.5]); hold on;
%set(gca,'Position',[0.4116    0.4583    0.1498    0.1425]);
set(gca,'Position',[0.4330 0.3900    0.1250    0.2150]);
xticks([-5 0 5]); xticklabels({'-5','0','5'});
yticks([-115 0 100 205]); yticklabels({'-115','0','100','205'});
colormap(gca,'default'); xlabel('x (cm)'); h=colorbar; ylabel(h,'$||m_{\mathrm{xy}}|/\mathrm{sin}(\alpha)-1|$', 'interpreter', 'latex', 'fontsize', 16);
plot(x(xmi:xma),b0map_x,'LineWidth',2.0,'Color',[0 0 0]); plot(x(xmi:xma),b0map_x+targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);
plot(x(xmi:xma),b0map_x-targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);plot(x(xmi:xma),b0map_x+max(fsim),'r--','LineWidth',1.0);
plot(x(xmi:xma),b0map_x+min(fsim),'r--','LineWidth',1.0); hold off; title(sprintf('Rel. Mag. Error at Y=%2.1f cm',y(y1)));  
% one y-location, phase error
subplot(336); im(x(xmi:xma),f,fliplr(squeeze(phaseE((nx-xma):(nx-xmi),y1,:))),[0 180]); hold on;
%set(gca,'Position',[0.6924    0.4583    0.1498    0.1425]);
set(gca,'Position',[0.7360 0.3900    0.1250    0.2150]);
xticks([-5 0 5]); xticklabels({'-5','0','5'});
yticks([-115 0 100 205]); yticklabels({'-115','0','100','205'});
colormap(gca,'default'); h=colorbar; ylabel(h,'$|<m_{xy}|(\circ)$','interpreter', 'latex', 'fontsize', 16);
plot(x(xmi:xma),b0map_x,'LineWidth',2.0,'Color',[0 0 0]); plot(x(xmi:xma),b0map_x+targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);
plot(x(xmi:xma),b0map_x-targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);plot(x(xmi:xma),b0map_x+max(fsim),'r--','LineWidth',1.0);
plot(x(xmi:xma),b0map_x+min(fsim),'r--','LineWidth',1.0); hold off; title(sprintf('Abs. Phase Error at Y=%2.1f cm',y(y1)));  

% one x-location, complex error
subplot(337); im(y(ymi:yma),f,fliplr(squeeze(complexE(x1,(ny-yma):(ny-ymi),:))),[0 3]); hold on;
%set(gca,'Position',[0.1308    0.1586    0.1429    0.1425]);
set(gca,'Position',[0.0700   0.0700    0.1850    0.2150]);
xticks([-7.5 0 10]); xticklabels({'-7.5','0','10'});
yticks([-115 0 100 205]); yticklabels({'-115','0','100','205'});
colormap(gca,'default'); ylabel('f (Hz)'); h=colorbar;  ylabel(h,'$|m_{\mathrm{xy}}/\mathrm{sin}(\alpha)-1|$', 'interpreter', 'latex', 'fontsize', 16);
plot(y(ymi:yma),b0map_y,'LineWidth',2.0,'Color',[0 0 0]); plot(y(ymi:yma),b0map_y+targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);
plot(y(ymi:yma),b0map_y-targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);plot(y(ymi:yma),b0map_y+max(fsim),'r--','LineWidth',1.0);
plot(y(ymi:yma),b0map_y+min(fsim),'r--','LineWidth',1.0,'Color',[1 0 0]); hold off; title(sprintf('Rel. Complex Error at X=%2.1f cm',x(x1)));  
% one x-location, magnitude error
subplot(338); im(y(ymi:yma),f,fliplr(squeeze(magE(x1,(ny-yma):(ny-ymi),:))),[0 1.5]); hold on;
%set(gca,'Position',[0.4116    0.1586    0.1429    0.1425]);
set(gca,'Position',[0.3730   0.0700    0.1850    0.2150]);
xticks([-7.5 0 10]); xticklabels({'-7.5','0','10'});
yticks([-115 0 100 205]); yticklabels({'-115','0','100','205'});
colormap(gca,'default'); xlabel('y (cm)'); h=colorbar;  ylabel(h,'$||m_{\mathrm{xy}}|/\mathrm{sin}(\alpha)-1|$', 'interpreter', 'latex', 'fontsize', 16);
plot(y(ymi:yma),b0map_y,'LineWidth',2.0,'Color',[0 0 0]); plot(y(ymi:yma),b0map_y+targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);
plot(y(ymi:yma),b0map_y-targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);plot(y(ymi:yma),b0map_y+max(fsim),'r--','LineWidth',1.0);
plot(y(ymi:yma), b0map_y+min(fsim),'r--','LineWidth',1.0); hold off; title(sprintf('Rel. Mag. Error at X=%2.1f cm',x(x1)));  
% one x-location, phase error
subplot(339); im(y(ymi:yma),f,fliplr(squeeze(phaseE(x1,(ny-yma):(ny-ymi),:))),[0 180]); hold on;
%set(gca,'Position',[0.6924    0.1586    0.1429    0.1425]);
set(gca,'Position',[0.6760   0.0700    0.1850    0.2150]);
xticks([-7.5 0 10]); xticklabels({'-7.5','0','10'});
yticks([-115 0 100 205]); yticklabels({'-115','0','100','205'});
colormap(gca,'default'); h=colorbar; ylabel(h,'$|<m_{xy}|(\circ)$','interpreter', 'latex', 'fontsize', 16);
plot(y(ymi:yma),b0map_y,'LineWidth',2.0,'Color',[0 0 0]); plot(y(ymi:yma),b0map_y+targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);
plot(y(ymi:yma),b0map_y-targetf/2,'--','LineWidth',1.0,'Color',[0.9 0.9 0.9]);plot(y(ymi:yma),b0map_y+max(fsim),'r--','LineWidth',1.0);
plot(y(ymi:yma),b0map_y+min(fsim),'r--','LineWidth',1.0,'Color',[1 0 0]); hold off; title(sprintf('Abs. Phase Error at X=%2.1f cm',x(x1)));  

return
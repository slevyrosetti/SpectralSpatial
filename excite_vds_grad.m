function [gx, gy, gz, kx, ky, kz, sx, sy, sz] = excite_vds_grad(ktype,FOV,gsamp,gmax,smax,npts,nil)
% [gx, gy, gz, kx, ky, kz, sx, sy, sz] = excite_vds_grad(ktype,FOV,gsamp,gmax,smax,npts,nil)
% Selects excitation gradient trajectory from various variable density
% spiral commands, includes gradient ramp-ups and rewinds
% By Sydney Williams, University of Michigan, 2016
%
% Inputs:
%   ktype       [1]             1: 1VDS at end 2: 1VDS at middle 3:1VDS at
%                               end 4: 2VDS low-res even 5: 2VDS low-res
%                               even alternating 6: 3VDS even 7: 3VDS even
%                               alternating 8: 2VDS med-res 9: 2VDS med-res
%                               alternating 10: 1VDS high-res
%   FOV         [1]             FOV_xy (cm)
%   gsamp       [1]             dwell time (ms)
%   gmax        [1]             max gradient amplitude (G/cm)
%   smax        [1]             max slew rate (G/cm/ms)
%   npts        [1]             number of time points in trajectory
%   nil         [1]             # of interleaves
%
% Outputs:
%   gx          [npts 1]        x gradient (G/cm)
%   gy          [npts 1]        y gradient (G/cm)
%   gz          [npts 1]        z gradient (G/cm)
%   kx          [npts 1]        x k-space (1/cm)
%   ky          [npts 1]        y k-space (1/cm)
%   kz          [npts 1]        z k-space (1/cm)
%   sx          [npts 1]        x slew (G/cm)
%   sy          [npts 1]        y slew (G/cm/ms)
%   sz          [npts 1]        z slew (G/cm/ms)
%
if ktype==1
    gshift=0;
    nk=1;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,153);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
            %g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    if gshift
        if gshift==1;
            g=circshift(g,[0 length(gout)]);
        else
            g=circshift(g,[0 (length(gout)+round(length(add0)/3))]);
        end
    else
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==2
    gshift=2;
    nk=1;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,153);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
            %g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    if gshift
        if gshift==1;
            g=circshift(g,[0 length(gout)]);
        else
            g=circshift(g,[0 (length(gout)+round(length(add0)/2))]);
        end
    else
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==3
    gshift=1;
    nk=1;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,153);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
            %g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    if gshift
        if gshift==1;
            g=circshift(g,[0 length(gout)]);
        else
            g=circshift(g,[0 (length(gout)+round(length(add0)/3))]);
        end
    else
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==4
    gshift=0;
    nk=2;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,153);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
            %g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    if gshift
        if gshift==1;
            g=circshift(g,[0 length(gout)]);
        else
            g=circshift(g,[0 (length(gout)+round(length(add0)/3))]);
        end
    else
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==5
    gshift=0;
    nk=2;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,153);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    if gshift
        if gshift==1;
            g=circshift(g,[0 length(gout)]);
        else
            g=circshift(g,[0 (length(gout)+round(length(add0)/3))]);
        end
    else
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==6
    gshift=0;
    nk=3;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,153);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
            %g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    if gshift
        if gshift==1;
            g=circshift(g,[0 length(gout)]);
        else
            g=circshift(g,[0 (length(gout)+round(length(add0)/3))]);
        end
    else
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==7
    gshift=0;
    nk=3;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,153);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    if gshift
        if gshift==1;
            g=circshift(g,[0 length(gout)]);
        else
            g=circshift(g,[0 (length(gout)+round(length(add0)/3))]);
        end
    else
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==8
    nk=2;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,300);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
            %g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms
elseif ktype==9
    nk=2;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,300);           % short >1ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    for i=1:nk
        if mod(i,2)==1
            g=[g add0 gout];            % adds zeros in center for each spiral, nk(i)
        elseif mod(i,2)==0
            g=[g add0 -gout];            % adds zeros in center for each spiral, nk(i) (rotating spiral)
        end
    end
    clear i;
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms 
elseif ktype==10
    nk=1;                                               % number of 2D k-space spirals
    xres=12.5;
    [ki,gi,si,ti] = vdsmex(nil,FOV,xres,gmax,smax*1e3,gsamp*1e-3,675);   % long 3ms spiral
    gi=[gi zeros(length(gi),1)];    % adds zero gz gradient
    gend = gi(end,:);               % rewind gradient
    ns = ceil(abs(gend)/smax/gsamp);% number of samples needed for rewind
    nsmax=max(ns);                  % max number of samples
    clear gdown;                    
    %rewinder
    for i=1:3                       % rewind down   
        gtemp = [ns(i)-1:-1:0]/ns(i)*gend(i);
        gtemp=padarray(gtemp,[0 abs(nsmax-ns(i))], 0, 'post');
        gdown(i,:)=gtemp;
    end
    gdown=gdown.';
    gout = [gi; gdown];             % appends rewind gradient
    ksum = sum(gout);               % units: grad-samples
    %triangle rephaser
    slpersamp = smax*gsamp;        % maximum slew rate per step down
    ns = ceil(sqrt(abs(ksum)/slpersamp)); % number of samples needed for rephaser
    nsmax=max(ns);                  % max number of samples
    clear tri;
  
    for i=1:3                       % triangle rephaser
        tritemp = [(1:ns(i)) (ns(i)-1:-1:0)];
        tritemp = tritemp.*ksum(i)/(sum(tritemp));
        tritemp=padarray(tritemp,[0 abs(2*nsmax-2*ns(i))], 0, 'post');
        tri(i,:)=tritemp;
    end
    gout=gout.';
    gout = [gout -tri];             % in G/cm
    rem=npts-nk*(length(gout)+1);   % remainder time samples left in trajectory
    add0=zeros(3,floor(rem/nk)-10); % determines number of 0s to add at k-space center
    gout=fliplr(gout);              % flips spiral so they are spiral in
    g=[];                           % blank array
    g=[g add0 gout];                % adds zeros in center for each spiral, nk(i)
    k = g2k(g',1)';
    s(1,:) = diff(g(1,:)./(gsamp));% slew rate vector in G/cm/ms
    s(2,:) = diff(g(2,:)./(gsamp));% slew rate vector in G/cm/ms
    s(3,:) = diff(g(3,:)./(gsamp));% slew rate vector in G/cm/ms
    tmax = length(g)*gsamp;       % in ms
    g=padarray(g,[0 (npts-length(g))], 0, 'post');
    k=padarray(k,[0 (npts-length(k))], 0, 'post');
    %k2=padarray(k2,[0 (npts-length(k))], 0, 'post');
    s=padarray(s,[0 (npts-length(s))], 0, 'post');
    kx=k(1,:)';                     % in 1/cm
    ky=k(2,:)';                     % in 1/cm
    kz=k(3,:)';
    gx=g(1,:)';                     % in G/cm
    gy=g(2,:)';                     % in G/cm
    gz=g(3,:)';                     % in G/cm
    sx=s(1,:)';                     % in G/cm/ms
    sy=s(2,:)';                     % in G/cm/ms
    sz=s(3,:)';                     % in G/cm/ms 
else                     % all gradient's off
    gx=zeros(npts,1);
    gy=zeros(npts,1);
    gz=zeros(npts,1);
    kx=zeros(npts,1);
    ky=zeros(npts,1);
    kz=zeros(npts,1);
    sx=zeros(npts,1);
    sy=zeros(npts,1);
    sz=zeros(npts,1);
end


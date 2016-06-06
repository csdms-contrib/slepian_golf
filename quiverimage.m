function quiverimage(data,lon,lat,del,renorm)
% quiverimage(data,lon,lat,del,renorm)
% 
% Plots the directions (quiver plot) that goes with IMAGEFNAN for
% tangential fields
%
% INPUT:
% 
% data      components in phi direction (data(:,:,1)) and in theta
%           direction (data(:,:,2))
% lon       longitude (phi) points for data
% lat       latitude (theta) points for data
% del       do not plot those vectors whose intensity is below del*maximum 
%           intensity value
% renorm    scale to maximum value (0) or multiply with the given value (if
%           not zero)
%
% See also IMAGEFNAN, QUIVPOLCAPS, QUIVERSPHERE
%
% Last modified by plattner-at-alumni.ethz.ch, 01/31/2012

defval('lon',linspace(0,360,size(data,2)));
defval('lat',linspace(90,-90,size(data,1)));
defval('del',0.01);
defval('renorm',0);

absdata=sqrt(data(:,:,1).^2 + data(:,:,2).^2);
maxab=max(max(absdata));
ind=find(absdata<del*maxab);
dat1=data(:,:,1); dat2=data(:,:,2);
[lon2,lat2]=meshgrid(lon,lat);
lon2(ind)=[]; lat2(ind)=[]; 
dat1(ind)=[]; dat2(ind)=[];
h=quiver(lon2,lat2,dat1,dat2,renorm,'ko'); 
set(h,'LineWidth',0.5)
set(h,'MarkerSize',1.5)
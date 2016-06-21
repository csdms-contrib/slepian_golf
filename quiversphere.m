function varargout=quiversphere(data,lon,lat,lift,del,renorm,lw,ms)
% [h,x,y,z,dx,dy,dz]=QUIVERSPHERE(data,lon,lat,lift,del,renorm,lw,ms)
%
% Plots tangential vector directions on globe. Use together with 
% PLOTONEARTH, as in the demos 
%
% INPUT:
% 
% data      The tangential vector field data to be plotted: 
%           A #lat x #lon x 2 matrix. data(:,:,1) are the lengths in phi-
%           direction and data(:,:,2) in theta-directions. For example use
%           the output of BLMCLM2XYZ
% lon       The coordinates of the longitude points 
% lat       The coordinates of the latitude points
%           Both, lon and lat must be in grid form (as for example output 
%           of matlab function MESHGRID)
% lift      Plot the strokes on the Sphere (0) or slightly above (value)
% del       Don't plot the strokes with magnitude smaller than del*max
%           amplitude.
% renorm    The 7th entry of the matlab QUIVER3 function. Set this to 0 if
%           the strokes should not be renormalized
% lw        Line width
% ms        Marker size
%
% OUTPUT:
%
% h         Graphic handle for the quiver plot
% x,y,z     3 dimensional coordinates of the stroke positions on the sphere
% dx,dy,dz  3 dimensional coordinates of the directions of the strokes at
%           each position on the sphere
%
% EXAMPLE: 
% 
% quiversphere('demo1') Plot the absolute values and the strokes for a 
%                       tangential Slepian for Eurasia on the globe
%
% quiversphere('demo2') Test if the tangential strokes are really
%                       tangential
%
% quiversphere('demo3') Plot tangential Slepian for Antarctica on the globe
% 
% quiversphere('demo4') Plot only on a part of the globe: Example Slepian 
%                       for Australia  
%
% See also PLOTONEARTH, PLOTONSPHERE, QUIVPOLCAPS, QUIVERIMAGE
%
% Last modified by plattner-at-alumni.ethz.ch 07/02/2013

defval('lon',[])
defval('lat',[])
defval('lift',0)
defval('del',0)
defval('renorm',1)
defval('lw',0.5)
defval('ms',1.5)

if ~isstr(data)

if ~isempty(lon)
  if size(data(:,:,1))~=size(lon) | size(data(:,:,1))~=size(lat)
    error('Wrong size arrays')
  end
end

nlon=size(data,2);
nlat=size(data,1);

if isempty(lon)
  % Make sphere for the data
  lons=linspace(0,360,nlon)/180*pi;
  lats=linspace(90,-90,nlat)/180*pi;
  [lons,lats]=meshgrid(lons,lats);
else
  lons=lon;
  lats=lat;
end

rads=ones(size(lats))+lift;
drads=zeros(size(data(:,:,1)));


[x,y,z]=sph2cart(lons,lats,rads);
%[dx,dy,dz]=vsph2vcart(lons,lats,rads,data(:,:,1),data(:,:,2),drads);
[dx,dy,dz]=dsph2dcart(lons,lats,data(:,:,1),data(:,:,2),drads);

% Remove vectors that are smaller than del*100% of the max abs vector
ab=sqrt(dx.^2+dy.^2+dz.^2);
maxab=max(max(ab));
ind=find(ab<del*maxab);
x(ind)=[]; y(ind)=[]; z(ind)=[]; dx(ind)=[]; dy(ind)=[]; dz(ind)=[];

h=quiver3(x,y,z,dx,dy,dz,renorm,'ok');
set(h,'LineWidth',lw)
set(h,'MarkerSize',ms)

axis image
view(140,30)

varns={h,x,y,z,dx,dy,dz};
varargout=varns(1:nargout);

elseif strcmp(data,'demo1')
    Lmax=24;
    index=100;% 15
    dom='eurasia'
    res=[1 3]%[0.1 2];
    comp='tangential';
    range=[0 360 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];
    [data,lat,lon,C,V,blmcosi,clmcosi]=vectorslepian(Lmax,...
        dom,comp,index,res,c11cmn);
    absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
    disp(sprintf('Eigenvalue %g',V(index)));
    dmax=max(max(absdata));
    absdata(find(absdata<dmax/100))=0;
    plotonearth(-absdata);
    caxis([-dmax dmax]);
    kelicol(1);
    hold on
    quiversphere(data{2},[],[],[],0.1);
    view(0,0)
    % Plot a circle around the world
    h=circ(1);
    % PlotonEarth has view(140,30). Need to also rotate the circle
    rotate(h,[1 0 0],90)
    rotate(h,[1 0 0],-30)
    rotate(h,[0 0 1],140)
    axis off
    % To make sure the circle is also rotated
    view(140,30)
    
    %fig2print(gcf,'fportrait')
    hold off
    %print('-depsc','-r600',sprintf('%s_%d_%d',dom,Lmax,index))
    
elseif strcmp(data,'demo2')
    % Test if the vectors are truly tangential
    Lmax=24;
    index=1;%15;
    dom='eurasia';%'namerica'
    res=[2 2];
    comp='tangential';
    range=[0 360 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];
    [data,lat,lon,C,V,blmcosi,clmcosi]=vectorslepian(Lmax,...
        dom,comp,index,res,c11cmn);
    [h,x,y,z,dx,dy,dz]=quiversphere(data{2}); 
    hold on
    inprod=x.*dx+y.*dy+z.*dz;
    plotonearth(inprod,1);
    disp(sprintf(...
        'Max inner product between tangential and radial vector is %g.',...
        max(max(abs(inprod)))));
    disp(sprintf(...
        'Tangential vector abs is %g',max(max(sqrt(dx.^2+dy.^2+dz.^2)))));    
    
    
elseif strcmp(data,'demo3')
    Lmax=18;
    index=11;
    dom='antarctica';%'namerica'%'eurasia'%
    res=[0.1 2];
    comp='tangential';
    range=[0 360 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];
    [data,lat,lon,C,V,blmcosi,clmcosi]=vectorslepian(Lmax,...
        dom,comp,index,res,c11cmn,[],[],1);
    absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
    dmax=max(max(absdata));
    absdata(find(absdata<dmax/100))=0;
    plotonearth(-absdata,1);
    caxis([-dmax dmax]);
    hold on
    quiversphere(data{2},[],[],[],0.01);
    kelicol(1)
    view(90,-90)
    axis tight
    hold off
    
    
elseif strcmp(data,'demo4')
    Lmax=20;
    index=3;
    dom='australia'
    res=[0.1 2];
    comp='tangential';
    off=10;
    XY=eval(sprintf('%s(10)',dom));
    range=[min(XY(:,1))-off max(XY(:,1))+off ...
        min(XY(:,2))-off max(XY(:,2))+off];
    c11cmn=[range(1) range(4) range(2) range(3)]; 
    [data,lat,lon,C,V,blmcosi,clmcosi]=vectorslepian(Lmax,...
        dom,comp,index,res,c11cmn);
    absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
    dmax=max(max(absdata));
    absdata(find(absdata<dmax/100))=0;
    [lons,lats]=meshgrid(lon{1},lat{1});
    plotonearth(-absdata,1,lons*pi/180,lats*pi/180);
    caxis([-dmax dmax]);
    hold on
    [lons2,lats2]=meshgrid(lon{2},lat{2});
    quiversphere(data{2},lons2*pi/180,lats2*pi/180,[],0.01);
    kelicol(1)
    axis tight
    hold off
             
end

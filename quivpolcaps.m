function varargout=quivpolcaps(data,del,renorm,lon,lat)
% [X,Y,Zx,Zy]=quivpolcaps(data,del,renorm,lon,lat)
%
% Quiver plot for polar cap regions transforms vectors using the
% transformation corresponding to PLOTPLM method 5 
%
% INPUT:
%
% data      data(:,:,1) is a 2d array containing the values of the vector 
%           field in longitudinal direction,(phi) 
%           data(:,:,2) is a 2d array containing the values of the vector 
%           field in latitudinal direction (theta)
% del       do not plot all vectors with length below del*(max abs value)
%           [default 0.05]
% renorm    scale to maximum value (0) or multiply with the given value (if
%           not zero)
% lon,lat   positions, in radians
%
% OUTPUT:
% 
% X     X Positions of the vectors in the XY coordinates looking from above
%       down to the polar cap
% Y     Y Positions of the vectors 
% Zx    Vector value in x-direction
% Zy    Vector value in y-direction
%
% See also DSPH2DCART, CAPVECTORSLEPIAN
%
% Last modified by plattner-at-alumni.ethz.ch, 05/09/2017

defval('lon',linspace(0,2*pi,size(data,2)));
defval('lat',linspace(pi/2,0,floor(size(data,1)/2)));% Only takes top half
defval('del',0.05)
defval('renorm',[]);
[LON,LAT]=meshgrid(lon,lat);

% Radius from 0 to 1; longitude is azimuth
r=cos(LAT);
x=r.*cos(LON);
y=r.*sin(LON);
Xpos=linspace(-1,1,length(lat));
Ypos=linspace(1,-1,length(lat));
y(1,:)=y(2,:)/2;
x(1,:)=x(2,:)/2;
phi=griddata(x,y,data(1:length(lat),:,1),Xpos,Ypos(:));
th =griddata(x,y,data(1:length(lat),:,2),Xpos,Ypos(:));

[X,Y]=meshgrid(Xpos,Ypos);

X=X(:); Y=Y(:); dphi=phi(:); dth=th(:); 
% Get PHI, TH for the new points
PHI=NaN(size(X));TH=NaN(size(Y));
for i=1:length(X)
    r=sqrt(X(i)^2 + Y(i)^2);
    TH(i)=acos(r);
    if(Y(i)>=0)
        PHI(i)=acos(X(i)/r);     
    else        
        PHI(i)=acos(-X(i)/r)+pi;
    end
end

[Zx,Zy,Zz]=dsph2dcart(PHI,TH,dphi,dth,zeros(size(PHI)));

% Here, take out the small values
absdata2=sqrt(Zx.^2 + Zy.^2);
maxab=max(max(absdata2));
ind=find(absdata2<del*maxab);
indNaN=find(isnan(absdata2));
ind=[ind;indNaN];
X(ind)=[]; Y(ind)=[]; Zx(ind)=[]; Zy(ind)=[];
if isempty(renorm)
    h=quiver(X,Y,Zx,Zy,'ko');
else
    h=quiver(X,Y,Zx,Zy,renorm,'ko');
end
xlim([-1 1]);
ylim([-1 1]);
axis square
set(h,'LineWidth',0.5)
set(h,'MarkerSize',1.5)

varns={X,Y,Zx,Zy};
varargout=varns(1:nargout);

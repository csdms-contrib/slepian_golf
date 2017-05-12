function varargout=GMTexportVec(Pcoef,Bcoef,Ccoef,fname,cmp,area,cutcirc,nanrng,res,onorout)
% [data,lon,lat]=GMTexportVec(Pcoef,Bcoef,Ccoef,fname,cmp,area,cutcirc,nanrng,res,onorout)
%
% Expands the vector spherical-harmonic field expressed by Pcoef, Bcoef,
% and Ccoef on a regular spatial grid and exports it in a GMT readable
% format. You can also only plot a limited area, or cut out a spherical
% cap.
%
% INPUT:
%
% Pcoef         The spherical harmonic coefficients for the Plm
% Bcoef         The spherical harmonic coefficients for the Blm
% Ccoef         The spherical harmonic coefficients for the Clm
% fname         Name for the output grid file
% cmpwrite      component: 
%               1 Radial (outward)
%               2 Southward (colatitudinal)
%               3 Eastward (longitudinal)
%               OR you can choose several: [1 2 3] (default: all)
% area          lon lat range for which to plot the field 
%               [lonmin latmax lonmax latmin]
% cutcirc       [TH clon ccola] parameters for the region you want to cut
%               out: TH is polar cap opening angle (onesided, like in 
%               Plattner & Simons (2014), ACHA), clon is longitude, ccola
%               is colatitude. ALL IN DEGREES
% nanrng        turn all values that are between -nanrng and nanrng into
%               nans (might help to make stuff invisible)
% res           Resolution in degrees per pixel
% onorout       Is the coef in ADDMON (0 or nothing), or ADDMOUT (1)?

%
% OUTPUT:
% 
% data          Calculated data, set to zero outside of cutcirc, and maybe
%               stuff turned into nans if nanrng is on
% lon,lat       Longitudes/latitudes of the model evaluation
%
% Last modified by plattner-at-alumni.ethz.ch, 5/9/2017

defval('cmp',[1 2 3])
defval('area',[0 90 360 -90])
defval('cutcirc',[])
defval('nanrng',0)
defval('res',0.2)
defval('onorout',0)

Lmax=sqrt(length(Pcoef))-1;

Plmcosi=coef2lmcosi(Pcoef,onorout);
Blmcosi=fcoef2flmcosi(Bcoef,onorout);
Clmcosi=fcoef2flmcosi(Ccoef,onorout);

if any(cmp==1)
    [data{1},lon,lat]=plm2xyz(Plmcosi,res,area);
else  
    data{1}=zeros(length(lat),length(lon));
end

if any(cmp==2) | any(cmp==3)
    [dat_t,lon,lat]=blmclm2xyz(Blmcosi,Clmcosi,res,area);  
    data{2}=dat_t(:,:,2);
    data{3}=dat_t(:,:,1);
    clear dat_t
else
    data{2}=zeros(length(lat),length(lon));
    data{3}=zeros(length(lat),length(lon));
end

rad=zeros(size(lon));

if ~isempty(cutcirc)
    cTH=cutcirc(1);
    clon=cutcirc(2);
    ccola=cutcirc(3);
    [data,lon,lat,rad,index]=cut2cap(data,lon,lat,rad,cTH,clon,clat);    
end

if nanrng   
   for cmpind=1:3
      indeggs=find(abs(data{cmpind})<=nanrng);
      data{cmpind}(indeggs)=nan;
   end
end

for cmpind=1:length(cmp)
    fprintf('Data range = %g to %g nT\n',...
        min(min(data{cmp(cmpind)})),max(max(data{cmp(cmpind)})))
    grdwrite2p(lon,flipud(lat),flipud(data{cmp(cmpind)}),...
        sprintf('%s_cmp%d.nc',fname,cmp(cmpind)));    
end


varns={data,lon,lat};
varargout=varns(1:nargout);

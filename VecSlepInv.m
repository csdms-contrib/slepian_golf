function varargout=VecSlepInv(data,theta,phi,dom,Lmax,Jrad,Jtan)
% [Pcoef,Bcoef,Ccoef,condi,dataweights]=VecSlepInv(data,theta,phi,dom,Lmax,Jrad,Jtan)
%
% Calculates best fitting vector spherical-harmonic coefficients Pcoef, 
% Bcoef, and Ccoef using iteratively reweighted least squares from the data
% within the region for maximum spherical-harmonic degree Lmax. 
%
% Future releases will allow for spherical caps rotated to other locations
% and for named regions. At this point, only polar caps
%
% INPUT:
%
% data      given as data{1}=rad component, data{2}=colat component,
%           data{3}=lon component, 
%           OR data is a vector of length 3*length(rad) as [rad;cola;lon]
% theta     colatitude of the data points in radians 0<=cola<=pi
% phi       longitude of the data points in radians 0<=phi<=2*pi
% dom       Semi-opening angle of the spherical cap or two angles for ring      
% Lmax      Maximum spherical harmonic degree
% Jrad      How many radial Slepian functions should be used to calculate
%           the solution? More means more sensitive to noise but higher 
%           spatial resolution
% Jtan      How many tangential Slepian functions should be used to calculate 
%           the solution? More means more sensitive to noise but higher
%           spatial resolution
%
% OUTPUT:
% 
%       !!!!ALL OUTPUT COEFFICIENTS ARE IN ADDMON FORMAT!!!
% Pcoef         Coefficients for the vector spherical-harmonics Plm
% Bcoef         Coefficients for the vector spherical-harmonics Blm
% Ccoef         Coefficients for the vector spherical-harmonics Clm
% condi         Conditioning number of the square matrix of evaluated 
%               Slepian functions
% dataweights   Final weights used in the iteratively reweighted least 
%               squares solution {1} radial, {2} tangential
%
% Last modified by plattner-at-alumni.ethz.ch, 5/10/2017

defval('niter',10)

if ~ischar(data)

%% Initialization
    
% Number of radial vector spherical-harmonics 
N_rad=(Lmax+1)^2;
% Number of tangential vector spherical-harmonics
N_tan=2*(Lmax+1)^2-2;
% Number of data points
npoints=length(theta);    
    
if length(data)==3
    datarad=data{1}(:);
    datatan=[data{2}(:);data{3}(:)];
    %data=[data{1};data{2};data{3}];
elseif length(data)==3*npoints
    datarad=data(1:npoints);
    datarad=datarad(:);
    datatan=data(npoints+1:end);
    datatan=datatan(:);
end

clear data;

%% Preparing the Slepian functions

% Get the coefficients. They don't need to be sorted. 
% We'll sort them after combining
[Gtan,Vtan]=vectanglmalpha(dom,Lmax,1);
if ischar(dom)
  [Grad,Vrad]=glmalpha(dom,Lmax);
else
  [Grad,Vrad]=glmalpharing(dom,Lmax,1);
end
% Only keen the first J
GradJ=Grad(:,1:Jrad);
GtanJ=Gtan(:,1:Jtan);
clear Grad; clear Gtan;

%% Solve the radial part 

% Calculate evaluated Plm functions. ylm is phase-shifted, 
% so we need the + pi
Plm=ylm([0 Lmax],[],theta,phi+pi,[],[],[],1);
% Calculate evaluated Slepian functions
M=GradJ'*Plm;

% Start with a least squares step:
MM=M*M';
Md=M*datarad;
slepcoef_rad=MM\Md;

% And now iteratively reweighted residual calculation
[slepcoef_rad,dataweights{1}]=itweighres(M,datarad,slepcoef_rad,niter);    

clear M;

% Now turn the resulting Slepian coefficients into vector
% spherical-harmonic coefficients
Pcoef=GradJ*slepcoef_rad;

%% Solve the tangential part

% Evaluate the vector spherical-harmonics
% The blmclm function is phase-shifted, so we need to include
[Blm,Clm]=blmclm(1:Lmax,[],theta,phi+pi,[],[],[],1);
% Calculate evaluated Slepian functions
M=GtanJ'*[Blm{1},Blm{2};Clm{1},Clm{2}];

% Start with a least squares step:
MM=M*M';
Md=M*datatan;
slepcoef_tan=MM\Md;

% And now iteratively reweighted residual calculation
[slepcoef_tan,dataweights{2}]=itweighres(M,datatan,slepcoef_tan,niter);    

clear M;

% Now turn the resulting Slepian coefficients into vector
% spherical-harmonic coefficients
Bcoef=GtanJ(1:N_tan/2,:)*slepcoef_tan;
Ccoef=GtanJ(N_tan/2+1:end,:)*slepcoef_tan;


%% Some finalization

% Now move them to addmon
Pcoef=out2on(Pcoef,Lmax);
Bcoef=out2on([0;Bcoef],Lmax);
Bcoef=Bcoef(2:end);
Ccoef=out2on([0;Ccoef],Lmax);
Ccoef=Ccoef(2:end);

if nargout>3
    condi=cond(MM);
else
    condi=[];
end

varns={Pcoef,Bcoef,Ccoef,condi,dataweights};
varargout=varns(1:nargout);


%% The demos
elseif strcmp(data,'demo1')

% Let's make random Pcoef, Bcoef, Ccoef, evaluate them, and then solve    
dom=[20,50];  % Semi opening angle of the polar cap
%dom='namerica'
Lmax=20;    % Maximum spherical-harmonic degree   
Jfactor=1;%2;  % How many times the Shannon number for inversion?
N=2000; % How many random data points?
noiselevel=0; % Percentage of noise?

res=0.5; % Plotting resolution

% Which component to show? 1=radial, 2=colatitudinal, 3=longitudinal
showcomp=3;

Jrad=round(Jfactor*spharea(max(dom),1)*(Lmax+1)^2);
Jtan=round(Jfactor*spharea(max(dom),1)*(2*(Lmax+1)^2-2));

% Make random model spherical harmonic coefficients or load them or
% simply write them down, or make random field following a power
% law. Here we just put in random for testing.
% This is one of Frederik's random field generator functions
Plmcosi_true=plm2rnd(Lmax,0);
Blmcosi_true=plm2rnd(Lmax,0);
% We will need to take away the L=0 part for the tangential vector
% spherical harmonics
Blmcosi_true=Blmcosi_true(2:end,:);
% Do the same for the Clm
Clmcosi_true=plm2rnd(Lmax,0);
Clmcosi_true=Clmcosi_true(2:end,:);
% Transform this into ADDMON coefficient vector
Pcoef_true=lmcosi2coef(Plmcosi_true,1);
Bcoef_true=flmcosi2fcoef(Blmcosi_true,1);
Ccoef_true=flmcosi2fcoef(Clmcosi_true,1);

% Make random points
[lon,lat]=randpatch(N,max(dom),0,0);
if length(dom)==2
    % if we want data in the ring, then just remove the points that are not
    % in the ring
    indskeep=(lat<90-min(dom));
    lon=lon(indskeep);
    lat=lat(indskeep);
end
    
% Evaluate the vector spherical harmonics
% This is for the radial
Plmeval=ylm([0 Lmax],[],(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1);

% And multiply with the coefficients to get the evaluated function
vecdata{1}=(Pcoef_true'*Plmeval)';
% Now for the tangential components
[Blmeval,Clmeval]=blmclm(1:Lmax,[],(90-lat)*pi/180,lon*pi/180+pi,[],[],[],1);
vecdata{2}=(Bcoef_true'*Blmeval{1} + Ccoef_true'*Clmeval{1})';
vecdata{3}=(Bcoef_true'*Blmeval{2} + Ccoef_true'*Clmeval{2})';

% Add noise 
for cmp=1:3
    vecdata{cmp}=vecdata{cmp} + noiselevel*mean(abs(vecdata{cmp}))*rand(size(vecdata{cmp})); 
end

% Now invert
[Pcoef_result,Bcoef_result,Ccoef_result]=VecSlepInv(vecdata,(90-lat)*pi/180,lon*pi/180,dom,Lmax,Jrad,Jtan);

% Now visually compare the results


%% visual comparison
% This is the data prep for plotting
switch showcomp
    case 1
        % The original data
        % We designed Pcoef_true as ADDMOUT, so it needs the 1
        Plmcosi_true=coef2lmcosi(Pcoef_true/sqrt(4*pi),1);
        [d_true,lonp,latp]=plm2xyz(Plmcosi_true,res);
        titlestring{1}='Radial original';
        
        % But the resulting coefficient is in ADDMON
        Plmcosi_result=coef2lmcosi(Pcoef_result/sqrt(4*pi));
        %plotplm(Plmcosi_result,[],[],2,res);
        [d_result,lonp,latp]=plm2xyz(Plmcosi_result,res);
        titlestring{2}='Radial result';
        
        % For comparison we have to turn Pcoef_true into addmon
        Pcoef_true_on=out2on(Pcoef_true,Lmax);
        Plmcosi_diff=coef2lmcosi(Pcoef_true_on-Pcoef_result);
        [d_diff,lonp,latp]=plm2xyz(Plmcosi_diff,res);
        titlestring{3}='Radial difference';
        
    case 2
        % The colatitudinal and longitudinal data require both Blm and Clm
        Blmcosi_true=fcoef2flmcosi(Bcoef_true/sqrt(4*pi),1);
        Clmcosi_true=fcoef2flmcosi(Ccoef_true/sqrt(4*pi),1);
        [dat,lonp,latp]=blmclm2xyz(Blmcosi_true,Clmcosi_true,res);
        % The colatitudinal data is the second part of the 3-D data matrix
        % See "help blmclm2xyz"
        d_true=dat(:,:,2);
        titlestring{1}='Colatitudinal original';
        
        % The colatitudinal and longitudinal data require both Blm and Clm
        Blmcosi_result=fcoef2flmcosi(Bcoef_result/sqrt(4*pi));
        Clmcosi_result=fcoef2flmcosi(Ccoef_result/sqrt(4*pi));
        [dat,lonp,latp]=blmclm2xyz(Blmcosi_result,Clmcosi_result,res);
        % The colatitudinal data is the second part of the 3-D data matrix
        % See "help blmclm2xyz"
        d_result=dat(:,:,2);
        titlestring{2}='Colatitudinal result';
        
        % For comparison we have to turn Bcoef_true and Ccoef_true into addmon
        % Remember the missing L=0. We need to add and then remove it to
        % use out2on
        Bcoef_true_on=out2on([0;Bcoef_true],Lmax);
        Ccoef_true_on=out2on([0;Ccoef_true],Lmax);
        Blmcosi_diff=fcoef2flmcosi(Bcoef_true_on(2:end)-Bcoef_result);
        Clmcosi_diff=fcoef2flmcosi(Ccoef_true_on(2:end)-Ccoef_result);        
        [dat,lonp,latp]=blmclm2xyz(Blmcosi_diff,Clmcosi_diff,res);
        % The colatitudinal data is the second part of the 3-D data matrix
        % See "help blmclm2xyz"
        d_diff=dat(:,:,2);
        titlestring{3}='Colatitudinal difference';
        
    case 3
        % Same as case 2 except we show other part of "3-D matrix"
        % The colatitudinal and longitudinal data require both Blm and Clm
        Blmcosi_true=fcoef2flmcosi(Bcoef_true/sqrt(4*pi),1);
        Clmcosi_true=fcoef2flmcosi(Ccoef_true/sqrt(4*pi),1);
        [dat,lonp,latp]=blmclm2xyz(Blmcosi_true,Clmcosi_true,res);
        % The colatitudinal data is the second part of the 3-D data matrix
        % See "help blmclm2xyz"
        d_true=dat(:,:,1);
        titlestring{1}='Longitudinal original';
        
        % The colatitudinal and longitudinal data require both Blm and Clm
        Blmcosi_result=fcoef2flmcosi(Bcoef_result/sqrt(4*pi));
        Clmcosi_result=fcoef2flmcosi(Ccoef_result/sqrt(4*pi));
        [dat,lonp,latp]=blmclm2xyz(Blmcosi_result,Clmcosi_result,res);
        % The colatitudinal data is the second part of the 3-D data matrix
        % See "help blmclm2xyz"
        d_result=dat(:,:,1);
        titlestring{2}='Longitudinal result';
        
        % For comparison we have to turn Bcoef_true and Ccoef_true into addmon
        % Remember the missing L=0. We need to add and then remove it to
        % use out2on
        Bcoef_true_on=out2on([0;Bcoef_true],Lmax);
        Ccoef_true_on=out2on([0;Ccoef_true],Lmax);
        Blmcosi_diff=fcoef2flmcosi(Bcoef_true_on(2:end)-Bcoef_result);
        Clmcosi_diff=fcoef2flmcosi(Ccoef_true_on(2:end)-Ccoef_result);
        [dat,lonp,latp]=blmclm2xyz(Blmcosi_diff,Clmcosi_diff,res);
        % The colatitudinal data is the second part of the 3-D data matrix
        % See "help blmclm2xyz"
        d_diff=dat(:,:,1);
        titlestring{3}='Longitudinal difference';
end
        
        
  
% Here is the actual plotting        

figure(1)
imagesc(lonp,latp,d_true);
axis xy
dmax=max(abs(caxis));
caxis([-1 1]*dmax)
kelicol(1)
colorbar
title(titlestring{1})

figure(2)
imagesc(lonp,latp,d_result);
axis xy
%dmax=max(abs(caxis));
caxis([-1 1]*dmax)
kelicol(1)
colorbar
title(titlestring{2})

figure(3)
imagesc(lonp,latp,d_diff);
axis xy
caxis([-1 1]*dmax)
kelicol(1)
colorbar
title(titlestring{3})

figure(4)
scatter(lon,lat,100,vecdata{showcomp},'.')
dmax=max(abs(caxis));
caxis([-1 1]*dmax)
xlim([0 360])
ylim([-90 90])
kelicol(1)
colorbar
title('original data')


end
    
    
    

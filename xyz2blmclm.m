function varargout=xyz2blmclm(fthph,L,method,lat,lon)
% [blmcosi,clmcosi]=xyz2blmclm(fthph,L,method,lat,lon)
%
% Forward vector spherical harmonic transform (decomposition) in the 4pi 
% normalized basis.
%
% Converts a spatially gridded field vector field into vector spherical 
% harmonics. For complete and regular spatial samplings [0 360 -90 90].
% If regularly spaced and complete, do not specify lat,lon.
% If not regularly spaced, fthph, lon and lat are column vectors.
%
% INPUT:
%
% fthph         tangential vector field defined on colatitude theta and 
%               longitude phi. fthph(:,:,1) for the phi direction component 
%               and fthph(:,:,2) for the theta direction component
% L             Maximum degree of the expansion (Nyquist checked)
% method        'gl'         By Gauss-Legendre integration (fast, accurate,
%                            preferred)
%               'simpson'    By Simpson integation (fast, inaccurate)
% lat           If not [90,-90], give latitudes explicitly, in degrees
% lon           If not [0,360], give longitudes explicitly, in degrees
%
% OUTPUT:
%
% blmcosi       Matrix listing l,m,cosine and sine coefficients for the blm
%               vector spherical harmonics
% clmcosi       Matrix listing l,m,cosine and sine coefficients for the clm
%               vector spherical harmonics
%
% Following the design of XYZ2PLM
%
% See also XYZ2PLM, VECTORSPECTRAL
%
%
% EXAMPLE:
% xyz2blmclm('demo1') Shows the performances for 'gl' and 'simpson', 
%                     randomly sampled in theta-direction for 'gl'                              
%
% Last modified by plattner-at-alumni.ethz.ch, 02/29/2012

defval('method','gl')
defval('lon',[])
defval('lat',[])
defval('dw',[])

if ~isstr(fthph)
    
as=0;
% If no grid is specified, assumes equal spacing and complete grid
if isempty(lat) & isempty(lon)
  % Decompose point vector values into phi and theta component
  fthph1=fthph(:,:,1);
  fthph2=fthph(:,:,2);
  clear fthph;
        
  % Test if data is 2D, and periodic over longitude
  fthph1=reduntest(fthph1);
  fthph2=reduntest(fthph2);
  
  polestest(fthph1);
  polestest(fthph2);
  % Make a complete grid
  nlon=size(fthph1,2);
  nlat=size(fthph1,1);
  % Nyquist wavelength
  Lnyq=min([ceil((nlon-1)/2) nlat-1]);
  % Colatitude and its increment
  theta=linspace(0,pi,nlat);
  as=1; % Equally spaced
  % Calculate latitude/longitude sampling interval; no wrap-around left
  dtheta=pi/(nlat-1);
  dphi=2*pi/nlon;
elseif isempty(lon)
  % If only latitudes are specified; make equal spacing longitude grid
  % Latitudes can be unequally spaced for 'gl'.
  
  % Decompose point vector values into phi and theta component
  fthph1=fthph(:,:,1);
  fthph2=fthph(:,:,2);
  clear fthph;
  fthph1=reduntest(fthph1);
  fthph2=reduntest(fthph2);
  theta=(90-lat)*pi/180;
  dtheta=(lat(1)-lat(2))*pi/180;
  nlat=length(lat);
  nlon=size(fthph1,2);
  dphi=2*pi/nlon;
  Lnyq=min([ceil((nlon-1)/2) ceil(pi/dtheta)]);
else
    error('Not yet implemented for irregularly sampled data')
end

% Decide on the Nyquist frequency
defval('L',Lnyq);
%disp(sprintf('Lnyq= %i ; expansion out to degree L= %i',Lnyq,L))

if L>Lnyq | nlat<(L+1)
  error('XYZ2BLMBLM: Function undersampled. Aliasing will occur.')
end

% Make cosine and sine matrices
[m,l,mz]=addmon(L);

blmcosi=[l m zeros(length(l),2)];
clmcosi=[l m zeros(length(l),2)];

% Define evaluation points
switch method
 case 'gl'
  % Highest degree of integrand will always be 2*L
  [w,x]=gausslegendrecof(2*L,[],[-1 1]);
  % Function interpolated at Gauss-Legendre latitudes; 2D no help
  fthph1=interp1(theta,fthph1,acos(x),'spline');
  fthph2=interp1(theta,fthph2,acos(x),'spline');  
 case {'simpson'}
  % Where else to evaluate the Legendre polynomials
  x=cos(theta);
 case {'irr','im'} 
  error('This method has not yet been implemented')        
 otherwise
  error('Specify valid method')
end 

fnpl=sprintf('%s/LSSM_TAN-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'LEGENDRE'),L,length(x));
 
if exist(fnpl,'file')==2 & as==1
 disp(sprintf('XYZ2BLMCLM Loading %s',fnpl))   
 load(fnpl)
else  
  % Evaluate Legendre polynomials and their derivatives at selected points

  if L>200
    h=waitbar(0,'Evaluating all Legendre polynomials');
  end   
 
  % Using the new ilk routine to calculate mPlm/sin and dPlm
  [mPlm,dPlm]=calcilk(L,x); 
  
  in1=0;
  in2=1;
  for l=0:L
    % Take the 1/sqrt(l*(l+1)) factor of the Blm Clm into account. Could do
    % that in calcilk
    mPlm(:,in1+1:in2)=1/sqrt(l*(l+1))*mPlm(:,in1+1:in2);
    dPlm(:,in1+1:in2)=1/sqrt(l*(l+1))*dPlm(:,in1+1:in2);   
    in1=in2;
    in2=in1+l+2;  
    if L>200
      waitbar((l+1)/(L+1),h)
    end
  end
  if L>200
    delete(h)
  end
  if as==1
    save(fnpl,'mPlm','dPlm')
  end
end
   
    % Perhaps demean the data for Fourier transform
    defval('dem',0)
    if strcmp(method,'im')&&bxon
      ep=0.8;
      badx=find(abs(sin(theta))<ep);
      ind=1:length(theta);
      ind=skip(ind,badx);
      fthph1=fthph1(ind,:);
      fthph2=fthph2(ind,:);
    end
    if dem
      meanm1=mean(fthph1,2);
      fthph1=fthph1-repmat(meanm1,1,nlon);
      meanm2=mean(fthph2,2);
      fthph2=fthph2-repmat(meanm2,1,nlon);
    end
    
    % Calculate integration over phi by the fast Fourier
    % transform. Integration of real input field with respect to the second
    % dimension of r, at  wavenumber m, thus at constant latitude. You get
    % as many wavenumbers m as there are longitudes; only use to L. With
    % Matlab's FFT, need to multiply by sampling interval.
    
    gfft1=dphi*fft(fthph1,nlon,2);
    gfft2=dphi*fft(fthph2,nlon,2);
    
    if dem
        % Add the azimuthal mean back in there
        gfft1(:,1)=2*pi*meanm1;
        gfft2(:,1)=2*pi*meanm2;
    end
    
    % Note these things are only half unique - the maximum m is nlon/2
    % But no Nyquist theory exists for the Legendre transform...
    cos_phi=real(gfft1);
    sin_phi=-imag(gfft1);
    cos_th=real(gfft2);
    sin_th=-imag(gfft2);    
    in1=0;
    in2=1;        

switch method    
 case 'gl'     
   in1=1;
   in2=3;
   %mPlm=Plm;
   % Loop over the degrees. Could go up to l=nlon if you want
   for l=1:L,
       % Integrate over theta using Gauss-Legendre integration        
       % Remember: Blm_th=dPlm, Clm_th=m*Plm, Blm_phi=m*Plm, Clm_phi=-dPlm
       % Be careful: for derivatives by phi, the sin/cos changes: sin-> cos
       % and cos ->-sin
       bphi_cos=sum(-sin_phi(:,1:l+1).*( diag(w)*mPlm(:,in1+1:in2)));
       bphi_sin=sum( cos_phi(:,1:l+1).*( diag(w)*mPlm(:,in1+1:in2)));
       cphi_cos=sum( cos_phi(:,1:l+1).*(-diag(w)*dPlm(:,in1+1:in2)));
       cphi_sin=sum( sin_phi(:,1:l+1).*(-diag(w)*dPlm(:,in1+1:in2)));
       
       bth_cos=sum( cos_th(:,1:l+1).*( diag(w)*dPlm(:,in1+1:in2)));
       bth_sin=sum( sin_th(:,1:l+1).*( diag(w)*dPlm(:,in1+1:in2)));
       cth_cos=sum(-sin_th(:,1:l+1).*( diag(w)*mPlm(:,in1+1:in2)));
       cth_sin=sum( cos_th(:,1:l+1).*( diag(w)*mPlm(:,in1+1:in2)));

       in1=in2;
       in2=in1+l+2;

       blmcosi(addmup(l-1)+1:addmup(l),3)=(bphi_cos(:)+bth_cos(:))/4/pi;
       blmcosi(addmup(l-1)+1:addmup(l),4)=(bphi_sin(:)+bth_sin(:))/4/pi;
       clmcosi(addmup(l-1)+1:addmup(l),3)=(cphi_cos(:)+cth_cos(:))/4/pi;
       clmcosi(addmup(l-1)+1:addmup(l),4)=(cphi_sin(:)+cth_sin(:))/4/pi;
   end
   rnk=[];
   blmcosi=blmcosi(2:end,:);
   clmcosi=clmcosi(2:end,:);
    
 case 'simpson'
  % Loop over the degrees. Could go up to l=nlon if you want
  for l=0:L,
    % Integrate over theta using Simpson's rule
    
    % Remember: Blm_th=dPlm, Clm_th=m*Plm, Blm_phi=m*Plm, Clm_phi=-dPlm
    % Be careful: for derivatives by phi, the sin/cos changes: sin-> cos
    % and cos ->-sin
    bphi_cos=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*(-sin_phi(:,1:l+1).*...
        mPlm(:,in1+1:in2)));
    bphi_sin=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*( cos_phi(:,1:l+1).*...
        mPlm(:,in1+1:in2)));     
    cphi_cos=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*(-cos_phi(:,1:l+1).*...
        dPlm(:,in1+1:in2)));  
    cphi_sin=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*(-sin_phi(:,1:l+1).*...
        dPlm(:,in1+1:in2)));  
    
    bth_cos=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*( cos_th(:,1:l+1).*...
        dPlm(:,in1+1:in2))); 
    bth_sin=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*( sin_th(:,1:l+1).*...
        dPlm(:,in1+1:in2))); 
    cth_cos=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*(-sin_th(:,1:l+1).*...
        mPlm(:,in1+1:in2))); 
    cth_sin=simpson(theta,...
        repmat(sin(theta(:)),1,l+1).*( cos_th(:,1:l+1).*...
        mPlm(:,in1+1:in2))); 

    in1=in2;
    in2=in1+l+2;
    % And stick it in a matrix [l m Ccos Csin]   
    blmcosi(addmup(l-1)+1:addmup(l),3)=(bphi_cos(:)+bth_cos(:))/4/pi;
    blmcosi(addmup(l-1)+1:addmup(l),4)=(bphi_sin(:)+bth_sin(:))/4/pi;
    clmcosi(addmup(l-1)+1:addmup(l),3)=(cphi_cos(:)+cth_cos(:))/4/pi;
    clmcosi(addmup(l-1)+1:addmup(l),4)=(cphi_sin(:)+cth_sin(:))/4/pi;
  end
  blmcosi=blmcosi(2:end,:);
  clmcosi=clmcosi(2:end,:);
    
    
end

% Get rid of machine precision error
blmcosi(abs(blmcosi(:,3))<eps,3)=0;
blmcosi(abs(blmcosi(:,4))<eps,4)=0;
clmcosi(abs(clmcosi(:,3))<eps,3)=0;
clmcosi(abs(clmcosi(:,4))<eps,4)=0;

varns={blmcosi,clmcosi};
varargout=varns(1:nargout);


elseif strcmp(fthph,'demo1')
    Lmax=5;
        
    [dems,dels,mz,lmc,mzin]=addmon(Lmax);
    cblm=rand((Lmax+1)^2,1);%zeros((Lmax+1)^2,1);%
    cclm=rand((Lmax+1)^2,1);%zeros((Lmax+1)^2,1);%
    cblm=reshape(insert(cblm,0,mzin),2,length(dems))';   
    cclm=reshape(insert(cclm,0,mzin),2,length(dems))';   
    dems=dems(2:end);
    dels=dels(2:end);
    cblm=cblm(2:end,:);
    cclm=cclm(2:end,:);    
    blmcosi=[dels dems cblm];
    clmcosi=[dels dems cclm];
    %blmcosi(2,4)=1;
    %blmcosi(8,4)=-1;
    %clmcosi(2,3)=-1;
     
    res=1;      
    range=[0 360 -90 90];
    c11cmn=[range(1) range(4) range(2) range(3)];

    [fthph,lon,lat]=blmclm2xyz(blmcosi,clmcosi,res); 
    % For the quiver plot, a coarser sampling
    [fthphb,lonb,latb]=blmclm2xyz(blmcosi,clmcosi,5); 

    % Plot the field we want to recover
    absfthph=sqrt(fthph(:,:,1).^2+fthph(:,:,2).^2);
    cax=[-max(max(absfthph)) max(max(absfthph))];
    imagefnan([0 90],[360 -90],absfthph,'kelicol',cax,[],1,100)     
    hold on   
    quiverimage(fthphb,lonb,latb)
    hold off
     
    % Set the number of sampling points for the GL recovery 
    Nlat=180;
    Nlon=360;
    % Set the random sampling points in latitude direction
    randpos=sort(rand(Nlat,1));
    lat=randpos*(range(3)-range(4))+range(4);
    lon=linspace(range(1),range(2),Nlon);
    [latg,long]=meshgrid(lat,lon);
    % Sample randomly in theta-direction, equidistantly in phi-direction
    fthphir=blmclm2xyz(blmcosi,clmcosi,latg(:),long(:));
    % Turn the sampled field values back into the right form
    ftpr=fthphir;
    clear fthphir;
    fthphir(:,:,1)=reshape(ftpr(:,:,1),Nlon,Nlat)';
    fthphir(:,:,2)=reshape(ftpr(:,:,2),Nlon,Nlat)';
    tic     
    % Recover th-randomly sampled field using GL
    [blmcosigl,clmcosigl]=xyz2blmclm(fthphir,Lmax,'gl',lat);
    time=toc;
    disp(sprintf('Calculation time for "gl" is %g sec',time))
    tic
    % Recover regularly sampled field using Simpson
    [blmcosisim,clmcosisim]=xyz2blmclm(fthph,Lmax,'simpson');
    time=toc;
    disp(sprintf('Calculation time for "simpson" is %g sec',time))
    xyz_sim=blmclm2xyz(blmcosisim,clmcosisim,1);
    xyz_gl=blmclm2xyz(blmcosigl,clmcosigl,1);
    
    % Plot everything
    figure
    cax=[-1 1];          
    subplot(2,3,1)
    imagefnan([],[],xyz_sim(:,:,1),[],cax)
    title('irr phi')
    subplot(2,3,2)
    imagefnan([],[],xyz_gl(:,:,1),[],cax)
    title('gl phi')
    subplot(2,3,3)
    imagefnan([],[],fthph(:,:,1),[],cax)
    title('true phi')
    subplot(2,3,4)
    imagefnan([],[],xyz_sim(:,:,2),[],cax)
    title('irr theta')
    subplot(2,3,5)
    imagefnan([],[],xyz_gl(:,:,2),[],cax)
    title('gl theta')
    subplot(2,3,6)
    imagefnan([],[],fthph(:,:,2),[],cax)
    title('true theta')
   
    figure
    subplot(2,3,1)
    imagesc(blmcosisim(:,3:4))
    title('Blm irr')
    colorbar
    subplot(2,3,2)
    imagesc(blmcosigl(:,3:4))
    title('Blm gl')
    colorbar
    subplot(2,3,3)
    imagesc(blmcosi(:,3:4))
    colorbar
    title('Blm true')
    subplot(2,3,4)
    imagesc(clmcosisim(:,3:4))
    title('Clm irr')
    colorbar
    subplot(2,3,5)
    imagesc(clmcosigl(:,3:4))
    title('Clm gl')
    colorbar
    subplot(2,3,6)
    imagesc(clmcosi(:,3:4))
    title('Clm true')
    colorbar
    
    disp(sprintf(...
   'Relative coefficient error of gl with random latitude points is %g',...
        sum(sum( (blmcosigl(:,3:4)-blmcosi(:,3:4)).^2 +...
        (clmcosigl(:,3:4)-clmcosi(:,3:4)).^2 ))/...
        sum(sum( (blmcosi(:,3:4)+clmcosi(:,3:4)).^2))     ));
    
     disp(sprintf(...
   'Relative coefficient error of simpson is %g',...
        sum(sum( (blmcosisim(:,3:4)-blmcosi(:,3:4)).^2 +...
        (clmcosisim(:,3:4)-clmcosi(:,3:4)).^2 ))/...
        sum(sum( (blmcosi(:,3:4)+clmcosi(:,3:4)).^2))     ));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grd=reduntest(grd)
% Tests if last longitude repeats last (0,360)
% and removes last data column
if sum(abs(grd(:,1)-grd(:,end))) >= size(grd,2)*eps*10
  disp(sprintf('Data violate wrap-around by %8.4e',...
		  sum(abs(grd(:,1)-grd(:,end)))))
end
grd=grd(:,1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function polestest(grd)
% Tests if poles (-90,90) are identical over longitudes 
var1=var(grd(1,:));
var2=var(grd(end,:));
if var1>eps*10 | var2>eps*10
  disp(sprintf('Poles violated by %8.4e and %8.4e',var1,var2))
end


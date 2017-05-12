function varargout=blmclm(l,m,theta,phi,check,tol,blox,irr)
% [B,C,theta,phi,dems,dels]=blmclm(l,m,theta,phi,check,tol,blox,irr)
%
% Calculates unit-normalized real vector spherical harmonics Blm, 
% DT (B.159). Also works for irregularly spaced lon/lat points.
%
% INPUT:
% 
% l      degree(s) (1 <= l <= infinity) [default: random]
% m      order (-l <= m <= l) [default: all orders -l:l]
%        l and m can be vectors, but not both at the sane time
%        if l is [0 L] and m is empty, do all degrees from 0 to L
%        if l is a single number and m is empty, it does m=0:l
% theta  colatitude vector (0<=theta<=pi) [default: 181 linearly spaced; 
%        not NaN!]
% phi    longitude vector (0<=phi<=2*pi) [default: 361 linearly spaced]
%        Unless irr=1, we assume you mean a 2-D (theta,phi) grid.
%        But if irr=1, length(theta(:)) must be equal to length(phi(:)).
% check  1 optional normalization check by Gauss-Legendre quadrature
%        0 no normalization check [default]
% tol    Tolerance for optional normalization checking [default: 1e-10]
% blox   0 Standard (lm) ordering, as ADDMOUT, l=0:L, m=-l:l [default]
%        1 Block-diagonal ordering, m=-L:L, l=abs(m):L
% irr    0 Regular grid, no matter how you input theta and phi [default]
%        1 Irregular grid, input interpreted as distinct pairs of theta, phi
%        Note that irr is a variable in the financial toolbox
%
% OUTPUT:
%
% B      The real vector spherical harmonics Blm at the desired argument(s):
%        B{1} (coefficients of the theta unit vector) and 
%        B{2} (coefficients of the phi unit vector).
%        Both as matrix with  dimensions of 
%           length(theta) x length(phi) x max(length(m),length(l)) 
%        OR 
%           (L+1)^2 x (length(theta)*length(phi)) if you put in
%           a degree l=[1:L] and an order []: lists orders -l to l 
%        OR
%            max(length(m),length(l)) or (L+1)^2 x length(theta) for irr=1
% C      The real vector spherical harmonics Clm at the desired argument(s),
%        exactly in the same format as the Blm
% theta  The latitude(s), which you might or not have specified
% phi    The longitude(s), which you might or not have specified
% dems   The orders to which the Blm and Clm belong
% dels   The degrees to which the Blm and Clm belong 
%        [if input needed interpreting]
%
% EXAMPLE:
%
% blmclm('demo1') Plots a tangential function using blmclm2xyz, using the
% coefficients times the values of blmclm, and compares the two
%
% See also LIBBRECHT, PLM, XLM, BLMCLM2XYZ, YLM
%
% Last modified by plattner-at-alumni.ethz.ch, 07/31/2012

% Default values
defval('l',round(rand*10))
defval('m',[])
defval('theta',linspace(0,pi,181))
defval('phi',linspace(0,2*pi,361))
defval('check',0)
defval('tol',-1)
defval('blox',0)
defval('irr',0)
dels=[];


if ~isstr(l)
  

  if blox~=0 & blox~=1
    error('Specify valid block-sorting option ''blox''')
  end

  if irr==1 & ~all(size(theta(:))==size(phi(:)))
    error('Input arrays must have the same dimensions for irregular grids')
  end


  % Make sure phi is a row vector
  phi=phi(:)';

  % If the degrees go from 1 to some L and m is empty, know what to do...
  if min(l)==1 & max(l)>1 & isempty(m)
    if irr==0
      % Here you assume a regular grid
      [PH,TH]=meshgrid(phi,theta);
    else
      PH=phi(:);
      TH=theta(:);
    end
    B{1}=repmat(NaN,(max(l)+1)^2,length(TH(:)));
    B{2}=repmat(NaN,(max(l)+1)^2,length(TH(:)));
    for thel=1:max(l)
      % Because here, too, you do the orders explicitly
      theB12=blmclm(thel,-thel:thel,theta,phi,check,tol,[],irr);      
      theB1=reshape(theB12{1},length(TH(:)),2*thel+1)';
      theB2=reshape(theB12{2},length(TH(:)),2*thel+1)';
      B{1}(thel^2+1:(thel+1)^2,:)=theB1;
      B{2}(thel^2+1:(thel+1)^2,:)=theB2;
    end  
    B{1}=B{1}(2:end,:);
    B{2}=B{2}(2:end,:);
    theta=TH(:);
    phi=PH(:);
    
    [dems,dels,mz,blkm]=addmout(max(l));
    if blox==1
      B{1}=B{1}(blkm,:);
      B{2}=B{2}(blkm,:);
      dems=dems(blkm);
      dels=dels(blkm);
    end
    C{1}= B{2};
    C{2}=-B{1};
    varns={B,C,theta,phi,dems,dels};
    varargout=varns(1:nargout);
    return
  end
  
  % ...if not, you're doing a single degree, perhaps as part of the above

  % Error handling common to PLM, XLM, YLM - note this resets defaults
  [l,m,mu,check,tol]=pxyerh(l,m,cos(theta),check,tol);

  % Straight to the calculation, check normalization on the XLM
  % Can make this faster, obviously, we're doing twice the work here
  [X,dX,theta,dems]=xdxlm(l,abs(m),theta,check);
  % The second dimension is always length(theta)
  
  
  if irr==0
    % Initialize the matrix with the spherical harmonics
    B{1}=repmat(NaN,[length(theta) length(phi) max(length(m),length(l))]);
    B{2}=repmat(NaN,[length(theta) length(phi) max(length(m),length(l))]);
    
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=repmat(diag(sqrt(2-(m(:)==0)))*...
	     cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))),length(l),1);
    dP=repmat(diag(sqrt(2-(m(:)==0)))*...
	   diag(m)*(-sin(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi)))),length(l),1);
    
    % Matrix of dimension max(length(l),length(m)) x length(theta)
    % Contains the 1./sin(theta) values and needs ti be element-wise
    % multiplied with X
    divsin=repmat(1./sin(theta(:)'),max(length(m),length(l)),1);
    
    % Matrix of dimension length(l) x length(theta)
    % Contains the 1/sqrt(l(l+1)) factor to be element-wise 
    % multiplied with X and dX
    % elle is a vector continaing either the different l or if we have a
    % single l and different m, then it is length(m) copies of the same l
    elle=repmat(l(:)',length(m),1);
    elle=elle(:)';
    divl=repmat(1./sqrt(elle.*(elle+1)),length(theta(:)),1)';
   
    
    % Make the real spherical harmonics
    % Here too, you assume a regular grid
    if prod(size(l))==1 & prod(size(m))==1
      B{1}=(divl.*dX)'*P;
      B{2}=(divl.*divsin.*X)'*dP;
      %Y=X'*P;
    else
      for index=1:max(length(m),length(l))
        B{1}(:,:,index)=(divl(index,:).*dX(index,:))'*P(index,:);
        B{2}(:,:,index)=(divl(index,:).*divsin(index,:).*X(index,:))'...
            *dP(index,:);
	    %Y(:,:,index)=X(index,:)'*P(index,:);
      end
    end
  else
   % Initialize the matrix with the spherical harmonics
    B{1}=repmat(NaN,[length(theta) max(length(m),length(l))]);
    B{2}=repmat(NaN,[length(theta) max(length(m),length(l))]);
    
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=repmat(diag(sqrt(2-(m(:)==0)))*...
	     cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))) ,length(l),1);
    % And also its derivative 
   dP=repmat(diag(sqrt(2-(m(:)==0)))*...
	   diag(m)*(-sin(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi)))),length(l),1);  
   
     % Matrix of dimension max(length(l),length(m)) x length(theta)
    % Contains the 1./sin(theta) values and needs ti be element-wise
    % multiplied with X
    divsin=repmat(1./sin(theta(:)'),max(length(m),length(l)),1);
    
    % Matrix of dimension length(l) x length(theta)
    % Contains the 1/sqrt(l(l+1)) factor to be element-wise 
    % multiplied with X and dX
    % elle is a vector continaing either the different l or if we have a
    % single l and different m, then it is length(m) copies of the same l
    elle=repmat(l(:)',length(m),1);
    elle=elle(:)';
    divl=repmat(1./sqrt(elle.*(elle+1)),length(theta(:)),1)';

    % Make the real spherical harmonics
    % Here too, you assume a regular grid
    if prod(size(l))==1 & prod(size(m))==1
      %Y=X.*P;
      B{1}=(divl.*dX).*P;
      B{2}=(divl.*divsin.*X).*dP;
    else
      for index=1:max(length(m),length(l))
        B{1}(:,index)=[(divl(index,:).*dX(index,:)).*P(index,:)]';
        B{2}(:,index)=[(divl(index,:).*divsin(index,:).*X(index,:))...
            .*dP(index,:)]';
	    %Y(:,index)=[X(index,:).*P(index,:)]';
      end
    end
  end
C{1}= B{2};
C{2}=-B{1};

varns={B,C,theta,phi,dems,dels};
varargout=varns(1:nargout);
  

% Examples

elseif strcmp(l,'demo1')
% Plot Blm and Clm using the classical way and then by evaluation
Lmax=18;
res=1;
reslow=5;
index=1;
dom='africa';

[B,C,theta,phi,dems,dels]=blmclm(Lmax);
[~,~,~,Cmat,Vmat,blmcosislp,clmcosislp]=vectorslepian(Lmax,dom,...
    'tangential',index);
coefs=blmclm2coef(blmcosislp,clmcosislp)';

% coefs=zeros(2*((Lmax+1)^2-1),1);
% coefs(7)=1;
% coefs([1 3 4]+end/2)=[0.3 0.2 0.5];

% First plot it the classical way
subplot(1,2,1)
[blmcosi,clmcosi]=coef2blmclm(coefs,Lmax);
[r,lon,lat]=blmclm2xyz(blmcosi,clmcosi,res,[0 89 360 -89]);
[rlow,lonplow,latplow]=blmclm2xyz(blmcosi,clmcosi,reslow,[0 89 360 -89]);
[lonlow,latlow]=meshgrid(lonplow,latplow);
absvals=sqrt(r(:,:,1).^2+r(:,:,2).^2);
plotplm(absvals,lon*pi/180,lat*pi/180,2);
clim=[-max(max(absvals)) max(max(absvals))];
caxis(clim);
hold on
quiversphere(rlow,lonlow,latlow);
hold off
kelicol(1);

% Now evaluate all Blm Clm at irregular points, multiply them with the 
% coeff and then plot
subplot(1,2,2)
% Use the same lat/lon as in the output of blmclm2xyz
lonp=lon;
latp=lat;
[lon,lat]=meshgrid(lonp,latp);
lat=90-lat;
la=lat*pi/180;
lo=lon*pi/180;
lo=lo(:);
la=la(:);

[B,C]=blmclm(1:Lmax,[],la,lo-pi,[],[],[],1);
% As in ylm, the output of blmclm is in the addmout format [-2 -1 0 1 2]. 
% But the coefficients are in the addmon format [0 -1 1 -2 2]. Luckilz, we
% have the out2on function, that rearanges the entries into the addmon.
% format

% To use out2on we need degrees 0 to Lmax. Easy: Simply add degree 0 with
% nans and then remove it again
Bon{1}=out2on([nan(1,size(B{1},2));B{1}],Lmax);
Bon{1}=Bon{1}(2:end,:);
Bon{2}=out2on([nan(1,size(B{2},2));B{2}],Lmax);
Bon{2}=Bon{2}(2:end,:);
Con{1}=out2on([nan(1,size(C{1},2));C{1}],Lmax);
Con{1}=Con{1}(2:end,:);
Con{2}=out2on([nan(1,size(C{2},2));C{2}],Lmax);
Con{2}=Con{2}(2:end,:);
fvals{1}=coefs'*([Bon{1};Con{1}])*sqrt(4*pi);
fvals{2}=coefs'*([Bon{2};Con{2}])*sqrt(4*pi);
newabsvals=sqrt(fvals{1}.^2+fvals{2}.^2);
data=reshape(newabsvals,length(latp),length(lonp));
plotplm(data,lonp*pi/180,latp*pi/180,2);
caxis(clim);

% Now the directions. Plot on wider grid
% lonplow=0:reslow:360;
% latplow=-89:reslow:89;
% [lonlow,latlow]=meshgrid(lonplow,latplow);
% Use the same lonlow latlow as in the output of blmclm2xyz
latlow=90-latlow;
lalow=latlow*pi/180;
lolow=lonlow*pi/180;
lolow=lolow(:);
lalow=lalow(:);

[Blow,Clow]=blmclm(1:Lmax,[],lalow,lolow-pi,[],[],[],1);
Bonlow{1}=out2on([nan(1,size(Blow{1},2));Blow{1}],Lmax);
Bonlow{1}=Bonlow{1}(2:end,:);
Bonlow{2}=out2on([nan(1,size(Blow{2},2));Blow{2}],Lmax);
Bonlow{2}=Bonlow{2}(2:end,:);
Conlow{1}=out2on([nan(1,size(Clow{1},2));Clow{1}],Lmax);
Conlow{1}=Conlow{1}(2:end,:);
Conlow{2}=out2on([nan(1,size(Clow{2},2));Clow{2}],Lmax);
Conlow{2}=Conlow{2}(2:end,:);
fvalslow{1}=coefs'*([Bonlow{1};Conlow{1}])*sqrt(4*pi);
fvalslow{2}=coefs'*([Bonlow{2};Conlow{2}])*sqrt(4*pi);
quivvals(:,:,1)=reshape(fvalslow{2},length(latplow),length(lonplow));
quivvals(:,:,2)=reshape(fvalslow{1},length(latplow),length(lonplow));
hold on
quiversphere(quivvals,lonlow,latlow);
hold off
kelicol(1);
max(max(max(abs(quivvals-rlow))))
absvals=reshape(absvals,1,size(absvals,1)*size(absvals,2));
max(abs(newabsvals-absvals))

end

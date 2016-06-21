function varargout=elm(Lmax,theta,phi,irr) 
% [E,theta,phi]=elm(Lmax,theta,phi,irr)
%
% Calculates unit-normalized gradient vector spherical harmonics Elm
% (see Handbook of Geomathematics, Plattner & Simons)
%
% INPUT:
% 
% Lmax      Maximum degree (we generate all Elm for L=0 to Lmax and all m)
% theta     colatitude for point evaluations 0<=theta<=pi
% phi       longitude for point evaluations 0<=phi<=2pi
% irr       0 for regular grid n=length(theta)*length(phi)
%           1 for irregular grid n=length(theta)=length(phi)
%
% OUTPUT:
% E         Gradient vector spherical harmonic evaluated at the n input 
%           points
%           E{1}(lm,n) radial component
%           E{2}(lm,n) theta component
%           E{3}(lm,n) phi component
% theta     evaluation point theta coordinates, length(theta)=n
% phi       evaluation points phi coordinates, length(phi)=n
%
% EXAMPLE: elm('demo1') to plot the rad and tan component of a single Elm
%
% Last modified by plattner-at-alumni.ethz.ch, 2/6/2013

defval('irr',0)

if ~isstr(Lmax)
if Lmax>1
    P=ylm([0 Lmax],[],theta,phi,[],[],[],irr);
    [B,~]=blmclm([1 Lmax],[],theta,phi,[],[],[],irr);
else
    error('Choose at least Lmax=2')
end


[~,~,~,~,~,~,~,bigl]=addmon(Lmax);


facP= sqrt((bigl+1)./(2*bigl+1));
facB=-sqrt( bigl(2:end)   ./(2*bigl(2:end)+1));
facPmat=spdiags(facP,0,length(bigl),length(bigl));
facBmat=spdiags(facB,0,length(bigl)-1,length(bigl)-1);

E{1}=facPmat*P;
E{2}=[zeros(1,size(P,2));facBmat*B{1}];
E{3}=[zeros(1,size(P,2));facBmat*B{2}];
% 
% 
% disp('WARNING: WRONG Elm!!! TESTING!!!')
% E{1}=facPmat*P;
% E{2}=[zeros(1,size(P,2));B{1}];
% E{3}=[zeros(1,size(P,2));B{2}];

varns={E,theta,phi};
varargout=varns(1:nargout);

elseif strcmp(Lmax,'demo1') 
   Lmax=5;
   ind=20;
   theta=0:0.05:pi;
   phi=0:0.05:2*pi;
   E=elm(Lmax,theta,phi,0);
   frad=reshape(E{1}(ind,:),length(theta),length(phi));
   fth=reshape(E{2}(ind,:),length(theta),length(phi));
   fph=reshape(E{3}(ind,:),length(theta),length(phi));
   subplot(3,1,1)
   imagesc(frad)
   subplot(3,1,2)
   imagesc(fth)
   subplot(3,1,3)
   imagesc(fph)
   
end
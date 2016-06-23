function varargout=pselm(TH,L)
% [lrnk,mrnk,lval,VV,Vsum]=PSELM(TH,L)
%
% Computes global ranking of eigenvalues for the concentration problem
% to the TANGENTIAL B-COMPONENT SINGLE SPHERICAL CAP considering all 
% angular orders m.
%
% INPUT:
%
% TH      Colatitudinal radius of the concentration region
% L       Bandwidth (maximum spherical harmonic degree)
% 
% OUTPUT:
% 
% lrnk    The index of the eigenvalue within its own m
% mrnk    The m of the eigenvalue
% lval    The sorted eigenvalues
% VV      The not globally sorted eigenvalue matrix (lrank-by-m)
% Vsum    Total sum of all the eigenvalues, including double counts
%
% Last modified by plattner-at-alumni.ethz.ch, 2/29/2012

defval('TH',40);
defval('L',18);

% Initialize eigenvalue matrix
VV=repmat(NaN,4*L,L+1);


% First m=0 because it has a differe
m=0;
Mmp=kerneltancapm(TH,L,m);
[~,Vp]=eig(Mmp);
VV(1:length(Vp),m+1)=diag(Vp);


for m=1:L;
  [Mmp,Mmm]=kerneltancapm(TH,L,m);
  [~,Vp]=eig(Mmp);
  [~,Vm]=eig(Mmm);
  V=[diag(Vp); diag(Vm)];
  VV(1:length(V),m+1)=V(:);
  % This is only filled for possible positive degrees from 0 to L
end

% Figure out GLOBAL rank ordering for the eigenvalues, with the repeated m
[a,b]=sort(VV(:));
% Put this in a matrix
mrnk=[repmat(0:(L+1),4*(L),1)];
lrnk=repmat([1:(L+1)]',1,4*(L));
b=b(~isnan(a)); b=flipud(b);
a=a(~isnan(a)); a=flipud(a);

% The ranked eigenvalues belong to this m
mrnk=mrnk(b);
% And they represent this number of nth eigenvalue
lrnk=lrnk(b);
lval=a;

% Figure out total sum of the eigenvalues including double counts
Vsum=VV;
Vsum(:,2:end)=Vsum(:,2:end);
Vsum=sum(Vsum(~isnan(Vsum)));

% Output
varn={lrnk,mrnk,lval,VV,Vsum};
varargout=varn(1:nargout);


function varargout=gradvecglmalphaptoJp(TH,L,phi,theta,omega,J)
% [Grot,V]=gradvecglmalphaptoJp(THL,phi,theta,omega,J)
%
% Gradient vectorial version if glmalphaptoJp. Loads glmalpha and rotates the first
% J Slepian functions only. 
%
% INPUT:
%
% TH        Angular extent of the spherical cap (degrees)
% L         Bandwidth (maximum angular degree) or passband (two degrees)
% phi       Longitude of the center (degrees) [default: 0]
% theta     Colatitude of the center (degrees) [default: 0]
% omega     Anticlockwise azimuthal rotation (degrees) [default: 0]
% J         The number of eigenfunctions that are being asked (and saved)
%
% OUTPUT:
%
% G        The unitary matrix of localization coefficients
% V        The eigenvalues 
%
% Last modified by plattner-at-alumni.ethz.ch, 02/11/2015

toout=1; % Should the result be in addmout?

defval('phi',0);
defval('theta',0);
defval('omega',0);

fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAPTOJP',...
         sprintf('gradvecglmalphaptoJp-%i-%i-%i-%i-%i-%i.mat',...
             TH,L,phi,theta,omega,J));

if exist(fname,'file')==2 
    load(fname)
    disp(sprintf('Loading %s',fname))
else
% Maybe this can be sped up by using the individual m blocks and not one
% big matrix
[H,S]=gradvecglmalphap(TH,L);
Hrot=nan((L+1)^2,J);
[S,isrt]=sort(S,2);
S=fliplr(S);
H=H(:,fliplr(isrt));
H=out2on(H(:,1:J),L);
[dems,dels,mz,lmcosi,mzi,mzo,bigm,bigl,rinm,ronm]=addmon(L);
%hh=waitbar(0,sprintf('Rotating the first %d Slepian functions',J));
disp('Calculating rotations in parallel mode')
try
    parpool
end
parfor j=1:J  
    lmcosip=lmcosi;  
    lmcosip(:,3:4)=reshape(insert(H(:,j),0,mzi),2,length(dems))';
    lmcosirot=plm2rot(lmcosip,omega,-theta,-phi);
    h=reshape(lmcosirot(:,3:4),1,2*length(dems));
    h=h(mzo);
    Hrot(:,j)=h(:); 
 %   waitbar(j/J,hh);
end
%delete(hh)
S=S(1:J);
% Transform back to addmout 
if toout
    Hrot=Hrot(rinm,:);
end
save(fname,'Hrot','S') 

end
varns={Hrot,S};
varargout=varns(1:nargout);
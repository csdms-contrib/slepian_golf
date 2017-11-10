function varargout=outgradvecglmalphatoJ(TH,L,phi,theta,omega,J)
% [Grot,V]=outgradvecglmalphatoJ(TH,L,phi,theta,omega,J)
%
% Outer field gradient vector Slepian function spherical cap rotation:
% Loads outgradvecglmalpha and rotates the first J Slepian functions only
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
% Grot     The unitary matrix of localization coefficients
% V        The eigenvalues 
%
% Last modified by plattner-at-alumni.ethz.ch, 02/25/2015

toout=1; % Should the result be in addmout?

defval('phi',0);
defval('theta',0);
defval('omega',0);

fname=fullfile(getenv('IFILES'),'OUTGRADVECGLMALPHATOJ',...
         sprintf('outgradvecglmalphauptoJp-%g-%i-%g-%g-%g-%i.mat',...
             TH,L,phi,theta,omega,J));

if exist(fname,'file')==2 
    load(fname)
    disp(sprintf('Loading %s',fname))
else         
    % Maybe this can be sped up by using the individual m blocks and not one
    % big matrix
    [H,S]=outgradvecglmalpha(TH,L);
    % H Should already be sparse but just to make sure:
    H=sparse(H);
    Hrot=sparse((L+1)^2-1,J);
    [S,isrt]=sort(S,2);
    S=fliplr(S);
    H=H(:,fliplr(isrt));
    
    %H=out2on(H(:,1:J),L);
    disp('Calculating rotations in parallel mode')
    try
        parpool
    end
    parfor j=1:J   
    %for j=1:J 
        lmcosip=[0 0 0 0;fcoef2flmcosi(H(:,j),1)];
        lmcosirot=plm2rot(lmcosip,omega,-theta,-phi);
        Hrot(:,j)=flmcosi2fcoef(lmcosirot(2:end,:),toout);
    end
    S=S(1:J);
    
    save(fname,'Hrot','S') 
    
end
varns={Hrot,S};
varargout=varns(1:nargout);    
    
    
    
         
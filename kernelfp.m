function KF=kernelfp(Lmax,dom,pars,ngl,rotb)
% KF=KERNELFP(Lmax,dom,pars,ngl,rotb)
%
% Calculates the vector Slepian kernel for the Flm. Can be calculated
% directly from the Kernels for the Plm (Ylm) and Blm which we already
% have.
%
% INPUT:
%
% Lmax      Maximum degree
% dom       Either the region 'england', 'eurasia', 'australia',
%           'greenland', 'africa', 'samerica', 'amazon', 'orinoco', 
%           'gpsnamerica', 'antarctica', 'alloceans', 'namerica' [default]
%           OR: 'patch' spherical patch with specs in 'pars'
%           OR: 'sqpatch' square patch with [thN thS phW phE] in 'pars'
%           OR: [lon lat] an ordered list defining a closed curve in
%           degrees
% pars      [th0,ph0,thR] for 'patch'
%                th0  Colatitude of the cap center, in radians
%                ph0  Longitude of the cap center, in radians
%                thR  Radius of the cap, in radians
%           OR: [thN thS phW phE] for 'sqpatch'
%           OR: N  splining smoothness for geographical regions 
%           [default: 10]
%           OR: the string with the name you want the result saved as
% ngl       an interger for Gauss-Legendre with the indicated number of
%           Gauss-Legendre points or 'paul' for Paul-Gaunt
% rotb      0 That's it, you got it [default: 0]
%           1 For, e.g., 'antarctica', if you were given rotated 
%           coordinates to make the integration procedure work, this option 
%           makes sure that the kernel matrix reflects this. If not, you 
%           have to apply counterrotation after diagonalizing in
%           LOCALIZATION.
%
% OUTPUT:
%
% KF        Localization kernel for the tangential plane
%           indexed as: degree  [1  1  1  2  2  2  2  2]
%                       order   [0 -1  1  0 -1  1 -2  2]
%           this is the ADDMON format
%
% See also KERNELEP, KERNELCP, KERNELBP, OUTGRADVECGLMALPHA
%
% Last modified by plattner-at-alumni.ethz.ch, 7/14/2017

defval('pars',[]);
defval('ngl',200)
defval('rotb',0)

try
    KP=kernelcp(Lmax,dom,pars,ngl,rotb);
    [~,KB]=kernelbp(Lmax,dom,pars,ngl,rotb);
catch
    KP=kernelc(Lmax,dom,pars,ngl,rotb);
    [~,KB]=kernelb(Lmax,dom,pars,ngl,rotb);
end

% Must remove degree 0 from KP:
KP=KP(2:end,2:end);

[~,~,~,~,~,~,~,bigl]=addmon(Lmax);

facB= sqrt((bigl(2:end)+1)./(2*bigl(2:end)+1));
facP= sqrt( bigl(2:end)   ./(2*bigl(2:end)+1));
facPmat=spdiags(facP,0,length(bigl)-1,length(bigl)-1);
facBmat=spdiags(facB,0,length(bigl)-1,length(bigl)-1);

KF_P=facPmat*KP*facPmat';
KF_B=facBmat*KB*facBmat';

KF=KF_P+KF_B;

% And remove numerical asymmetries
KF=(KF+KF')/2;
KF=full(KF);



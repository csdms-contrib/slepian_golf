function KT=kerneltorp(Lmax,dom,pars,ngl,rotb)
%  KT=kerneltorp(Lmax,dom,pars,ngl,rotb)
%
% Calculates the vector Slepian kernel for the Clm (torroidal part of the 
% field). Can be calculated directly from the Kernels of Blm which we already
% have.
%
% It is exactly the same as the kernel for the Blm, soThis function is
% simply a wrapper function for kernelb or kernelbp.
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
% K         Localization kernel for the tangential plane
%           indexed as: degree  [1  1  1  2  2  2  2  2]
%                       order   [0 -1  1  0 -1  1 -2  2]
%           this is the ADDMON format
%
% See also KERNELFP, KERNELEP, KERNELCP, KERNELBP, OUTGRADVECGLMALPHA
%
% Last modified by plattner-at-alumni.ethz.ch, 2/25/2015

defval('pars',[]);
defval('ngl',200)
defval('rotb',0)

try
    [~,KT]=kernelbp(Lmax,dom,pars,ngl,rotb);
catch
    [~,KT]=kernelb(Lmax,dom,pars,ngl,rotb);
end





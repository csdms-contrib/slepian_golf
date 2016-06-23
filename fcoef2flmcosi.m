function flmcosi=fcoef2flmcosi(fcoef,onorout)
% flmcosi=fcoef2flmcosi(fcoef,onorout)
%
% Transform coefficients for outer vector spherical harmonics Flm or
% toroidal spherical harmonics Clm
% given as a vector "fcoef" into the flmcosi format.
% Requires the fcoef vector to be "complete" as in having (L+1)^2-1 entries.
%
% INPUT:
%
% fcoef     the spherical-harmonic coefficients for Flm or Clm ordered in 
%           either addmon or addmout format.
% onorout   if fcoef in addmout format, set this to 1. Otherwise (addmon) 
%           leave out or set it to zero
% 
% OUTPUT:
%
% flmcosi    the coefficients ordered in lmcosi format (l m cos sin)
%
% See also  coef2lmcosi, lmcosi2coef, coef2blmclm, blmclm2coef
%
% Last modified by plattner-at-alumni.ethz.ch, 2/23/2015

defval('onorout',0)

fcoef=[0;fcoef(:)];
flmcosi=coef2lmcosi(fcoef,onorout);
flmcosi=flmcosi(2:end,:);

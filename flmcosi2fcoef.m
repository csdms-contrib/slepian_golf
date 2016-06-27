function fcoef=flmcosi2fcoef(flmcosi,onorout)
% fcoef=flmcosi2fcoef(flmcosi,onorout)
%
% Transform coefficients for outer vector spherical harmonics Flm or the 
% toroidal vector spherical harmonics Cl mgiven in flmcosi format 
% (= lmcosi without the L=0 part) into
% vector "fcoef" (without L=0) either in ADDMON or ADDMOUT.
%
% INPUT:
%
% flmcosi   the coefficients ordered in lmcosi format (l m cos sin)
% onorout   if you want fcoef in addmout format, set this to 1. 
%           Otherwise (addmon) leave out or set it to zero
%
% OUTPUT:
%
% fcoef     the spherical-harmonic coefficients for Flm or Clm ordered in 
%           either addmon or addmout format.
%
% See also  fcoef2flmcosi, coef2lmcosi, lmcosi2coef, 
%           coef2blmclm, blmclm2coef
%
% Last modified by plattner-at-alumni.ethz.ch, 6/27/2016

defval('onorout',0)

fcoef=lmcosi2coef([0 0 0 0;flmcosi],onorout);

fcoef=fcoef(2:end);

function coef=blmclm2coef(blmcosi,clmcosi)
% coef=blmclm2coef(blmcosi,clmcosi)
%
% Transforms blmcosi and clmcosi into a vector of coefficients that can be
% fed into vecslepanalysis
%
% INPUT:
%
% blmcosi   The coefficients for the Blm vector spherical harmonics, 
%           prepared for BLMCLM2XYZ
% clmcosi   The coefficients for the Clm vector spherical harmonics, 
%           prepared for BLMCLM2XYZ
%
% OUTPUT: 
%
% coef      Vector spherical harmonics coefficient vector, that can be
%           used for VECSLEPANALYSIS 
%
% Last modified by plattner-at-alumni.ethz.ch, 02/29/2012
%
% See also COEF2BLMCLM

Lmax=blmcosi(end,1);

[demsz,delsz,mz,lmc,mzin,mzo]=addmon(Lmax); 

% Add the l=0 part to make it easier
blmcosi=[NaN(1,4);blmcosi];
blm=reshape(blmcosi(:,3:4),1,2*length(demsz));
bcoeff=blm(mzo);
bcoeff=bcoeff(2:end);
if ~isempty(clmcosi)
    clmcosi=[NaN(1,4);clmcosi];
    clm=reshape(clmcosi(:,3:4),1,2*length(demsz));
    ccoeff=clm(mzo);
    ccoeff=ccoeff(2:end);
else
    ccoeff=[];
end
    
coef=[bcoeff ccoeff];


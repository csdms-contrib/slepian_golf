function [blmcosi,clmcosi]=coef2blmclm(coef,Lmax)
% [blmcosi,clmcosi]=COEF2BLMCLM(coef,Lmax)
% 
% Turns a coefficient vector, as for example an eigenvector of a tangential
% Slepian concentration kernel matrix from KERNELBP, into blmcosi and 
% clmcosi, as they are required for BLMCLM2XYZ. 
%
% INPUT:
% 
% coef      Vector spherical harmonics coefficient vector, as for example
%           the eigenvector values, sorted as in the matrix that comes out
%           of KERNELB or KERNELBP
% Lmax      The maximum degree
%
% OUTPUT:
%
% blmcosi   The coefficients for the Blm vector spherical harmonics, 
%           prepared for BLMCLM2XYZ
% clmcosi   The coefficients for the Clm vector spherical harmonics, 
%           prepared for BLMCLM2XYZ
%
% Last modified by plattner-at-alumni.ethz.ch, 08/2/2012
%
% See also KERNELB, KERNELBP, BLMCLM2XYZ, BLMCLM2COEF

[demsz,delsz,mz,lmc,mzin]=addmon(Lmax);
coef=coef(:);

% First, only take those coefficients that belong to the blm
coefblm=coef(1:size(coef,1)/2);
% And those that only belong to the clm
coefclm=coef((size(coef,1)/2+1):end);    
% The tangential field does not have an l=0 coefficient. In order to
% still use the already available functions, proceed as follows
% First put in dummy values for l=0 and then remove the l=0 part again
dems=demsz(2:end);
dels=delsz(2:end);

coefblm=[NaN;coefblm];    
coefblm=reshape(insert(coefblm,0,mzin),2,length(demsz))';  
coefblm=coefblm(2:end,:);    
blmcosi=[dels dems coefblm];

coefclm=[NaN;coefclm];    
coefclm=reshape(insert(coefclm,0,mzin),2,length(demsz))';   
coefclm=coefclm(2:end,:); 
clmcosi=[dels dems coefclm];


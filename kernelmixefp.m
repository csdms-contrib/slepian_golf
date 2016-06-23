function KEF=kernelmixefp(Lin,Lout,dom,pars,ngl,rotb)
% KEF=kernelmixefp(Lin,Lout,dom,pars,ngl,rotb)
%
% Calculates the mixed product E*F part of the overal inner+outer altitude
% vector Slepian kernel. Reuses the kernels for the Plm (Ylm) and Blm which
% we already have.
%
% INPUT:
%
% Lin       Maximum degree for the inner source component of the altitude
%           vector Slepian functions
% Lout      Maximum degree for the outer source component of the altitude
%           vector Slepian functions
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
% KEF       Localization kernel for the tangential plane
%           indexed as: degree  [1  1  1  2  2  2  2  2]
%                       order   [0 -1  1  0 -1  1 -2  2]
%           this is the ADDMON format
%
% See also KERNELEP, KERNELFP, KERNELCP, KERNELBP, OUTGRADVECGLMALPHA
%
% Last modified by plattner-at-alumni.ethz.ch, 4/15/2015

defval('pars',[]);
defval('ngl',200)
defval('rotb',0)

% We will load the entire kernels for the larger max degree because if you
% are doing this you will need to calculate the entire kernels anyhow for
% the ElmElm or FlmFlm part
Lmax=max(Lin,Lout);

try
    KP=kernelcp(Lmax,dom,pars,ngl,rotb);
    [~,KB]=kernelbp(Lmax,dom,pars,ngl,rotb);
catch
    KP=kernelc(Lmax,dom,pars,ngl,rotb);
    [~,KB]=kernelb(Lmax,dom,pars,ngl,rotb);
end

% First make sure that the right ls are covered. Remember: The columns
% represent the ls for Flm, so no l=0 and the rows are for the coefficients
% for the Elm, so we need the l=0.
% So: Remove the first column of KP and add a zero row on top of KB

KP=KP(:,2:end);
KB=[zeros(1,size(KB,2));KB];
% They should now have same size.

% Now only take the parts that we really need. Remember, here we don't know
% if Lin or Lout is greater
KP=KP(1:(Lin+1)^2,1:(Lout+1)^2-1);
KB=KB(1:(Lin+1)^2,1:(Lout+1)^2-1);

[~,~,~,~,~,~,~,biglE]=addmon(Lin);
[~,~,~,~,~,~,~,biglF]=addmon(Lout);

facPE= sqrt((biglE(1:end)+1)./(2*biglE(1:end)+1));
facPF= sqrt( biglF(2:end)   ./(2*biglF(2:end)+1));
facBE=-sqrt( biglE(1:end)   ./(2*biglE(1:end)+1));
facBF= sqrt((biglF(2:end)+1)./(2*biglF(2:end)+1));

facPEmat=spdiags(facPE,0,length(biglE),length(biglE));
facPFmat=spdiags(facPF,0,length(biglF)-1,length(biglF)-1);
facBEmat=spdiags(facBE,0,length(biglE),length(biglE));
facBFmat=spdiags(facBF,0,length(biglF)-1,length(biglF)-1);

KEF_P=facPEmat*KP*facPFmat';
KEF_B=facBEmat*KB*facBFmat';

KEF=KEF_P+KEF_B;

KEF=full(KEF);



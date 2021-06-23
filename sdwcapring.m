function [Vm,Cm]=sdwcapring(TH,L,m)
% [Vm,Cm]=sdwcapring(TH,L,m)
%
% Returns the scalar/radial Slepian coefficients for a single order m block
% for spherical caps or rings.
%
% INPUT:
%
% TH    Semi-opening angle of the spherical cap or two semi-opening angles
%       for the ring between them
% L     Maximum spherical-harmonic degree
% m     spherical-harmonic order
% 
% OUTPUT:
%
% Gm    Matrix, who's columns are the eigenvectors = spherical-harmonic
%       coefficients for the Slepian functions for single order m
% Vm    Vector of eigenvalues in the same order as the columns of G for
%       single order m
%
% Last modified by plattner-at-alumni.ethz.ch, 5/11/2017


% See if already calculated
dirname=fullfile(getenv('IFILES'),'SDWCAPRING');
fnpl=fullfile(dirname,sprintf(...
      'SDWRING-%f-%f-%i-%i.mat',max(TH),min(TH),L,m));
if exist(fnpl,'file')==2
    load(fnpl)
    disp(sprintf('%s loaded',fnpl))
else  

    % Set up K
    if length(TH)==2 
        K1=kernelcapring(max(TH),L,m);
        K2=kernelcapring(min(TH),L,m);        
        K=K1-K2;        
    elseif length(TH)==1
        K=kernelcapring(TH,L,m);        
    end

    % Now eigendecomposition
    [Cm,Vm]=eig(K);
    Vm=diag(Vm);

    try
    	% Matlab
    	save(fnpl,'Cm','Vm','-v7.3')
    catch
    	% Octave
    	save(fnpl,'Cm','Vm')      
    end

end % End of calculation

end % end of function
    




% Function to set up single K
function K=kernelcapring(TH,L,m)
% This I took from sdwcap.m
% Convert to radians
TH=TH*pi/180;
maxL=max(L);

% Setting Bandpass to zero for now
bp=0;
% Initialize kernel
K=zeros(maxL+1-max(m,bp*min(L)),maxL+1-max(m,bp*min(L)));

% Construct kernel: integral of Ylm over the cap
lmin=max(m,bp*min(L));

for lr=lmin:maxL
  for lc=lr:maxL
    % Orthonormalization of Ylm is to unity over unit sphere
    % When TH=180, K should be the identity matrix
    % The pi*(1+(m==0)) is the longitudinal integral
    % Note that this line ALSO would work with 'paul' but then it
    % would be very inefficient in not reusing the database
    K(lr+1-lmin,lc+1-lmin)=...
        legendreprodint(lr,m,lc,m,cos(TH),'gl')...
        *sqrt(2*lr+1)*sqrt(2*lc+1)/(4*pi)*pi*(1+(m==0));
    % Symmetrize
    K(lc+1-lmin,lr+1-lmin)=...
        K(lr+1-lmin,lc+1-lmin);
  end
end

end

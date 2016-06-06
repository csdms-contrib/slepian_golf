function [P,B]=EFC2PBC(E,F,Lmax)
% [P,B]=EFC2PBC(E,F,Lmax)
%
% Transforms the vector spherical harmonic coefficients from the EFC basis
% (up/downward continuable) to the PBC basis (classical). 
% 
% Two different modes: 
% EITHER the E and F coefficients are given as Elmcosi, Flmcosi;
% then the returning P B will also be in the lmcosi format (Plmcosi, 
% Blmcosi),
% OR E and F are given in either the ADDMOUT or ADDMON format 
% (both work because the factor depends on the degree, not on the order).
% In this case the returning coefficients will be given in the same format
%
% Remarks: 1) C is not required, because it is not transformed
%          2) E and P have degrees 0,1,..., F and B have degrees 1,... 
%
% INPUT: 
%
% E, F  up/downward continuable vector spherical harmonics coefficients. 
%       Either in the lmcosi or in the ADDMOUT or ADDMON format. If given 
%       in the ADDMOUT or ADDMON format, Lmax is required and all 
%       coefficients are needed.
% Lmax  Maximum degree. Only needed if the coefficients are given in either
%       the ADDMON or ADDMOUT format. Must not be provided if lmcosi format
%       is used.
%
% OUTPUT:
%
% P, B  Classical vector spherical harmonics coefficients (radial and 
%       projected gradient tangential).
%       Either given as lmcosi or vector of coefficients (depending on the 
%       input).
%
% EXAMPLE:
%
% EFC2PBC('demo1')  Generate random P B, transform to E F and back
%
% COMMENT:
%
% This is done using a for loop. I guess it would be faster if the loop is
% replaced by a vector multiplication.
%
% Last modified by plattner-at-alumni.ethz.ch, 08/20/2012
%
% See also PBC2EFC


defval('Lmax',-1);

if ~ischar(E)

% Decision if P and B are given as lmcosi
if (size(E,2)==4)&&(size(F,2)==4)&&(Lmax==-1)
    % The lmcosi case 
    
    % Remove the E_00 coefficient from the list to make things easier
    Eoo=E(1,:);
    E=E(2:end,:);
    
    P=nan(size(E)); B=nan(size(F));
    % Copy the l and m values
    P(:,1)=E(:,1); P(:,2)=E(:,2);
    B(:,1)=F(:,1); B(:,2)=F(:,2);
    
    for i=1:size(E,1)
        L=E(i,1);
        Pfac1= sqrt((L+1)/(2*L+1));
        Pfac2= sqrt((L  )/(2*L+1));
        Bfac1=-Pfac2;
        Bfac2= Pfac1;

        P(i,3)=Pfac1*E(i,3) + Pfac2*F(i,3); 
        P(i,4)=Pfac1*E(i,4) + Pfac2*F(i,4);

        B(i,3)=Bfac1*E(i,3) + Bfac2*F(i,3); 
        B(i,4)=Bfac1*E(i,4) + Bfac2*F(i,4); 

    end
    
    % And now put the P_00 coefficient into E
    P=[Eoo;P];
    
elseif (size(E,1)==(Lmax+1)^2) && (size(F,1)==(Lmax+1)^2-1)
    % The list case    
    
    % Remove the E_00 coefficient to make things easier
    Eoo=E(1,:);
    E=E(2:end,:);
    
    [~,~,~,~,~,~,~,bigl]=addmon(Lmax);
    bigl=bigl(2:end);
    P=nan(size(E)); B=nan(size(F));
    for i=1:(Lmax+1)^2-1
        L=bigl(i);
        Pfac1= sqrt((L+1)/(2*L+1));
        Pfac2= sqrt((L  )/(2*L+1));
        Bfac1=-Pfac2;
        Bfac2= Pfac1;

        P(i,:)=Pfac1*E(i,:) + Pfac2*F(i,:); 
        B(i,:)=Bfac1*E(i,:) + Bfac2*F(i,:);
    end
    
    % And now put the P_00 coefficient into E
    P=[Eoo;P];
    
else
    error('Something is wrong with the input. Or you use vector and forgot to provide Lmax')
end    
  


elseif strcmp(E,'demo1')
% Start with P B, transform to E F and back for both, lmcosi and vector    
Lmax=5;
Plmcosi=plm2rnd(Lmax);

Blmcosi=plm2rnd(Lmax);
Blmcosi=Blmcosi(2:end,:);    

vec=blmclm2coef(Plmcosi(2:end,:),Blmcosi);
vec=vec(:);
P=[Plmcosi(1,3);vec(1:end/2)]; B=vec(end/2+1:end);

% Elmcosi, Flmcosi
[Elmcosi,Flmcosi]=PBC2EFC(Plmcosi,Blmcosi);
[Plmcosi2,Blmcosi2]=EFC2PBC(Elmcosi,Flmcosi);

fprintf('Difference after double transformation %g\n',...
    norm(Plmcosi(:)-Plmcosi2(:))+norm(Blmcosi(:)-Blmcosi2(:)));


end
    

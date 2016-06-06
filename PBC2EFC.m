function [E,F]=PBC2EFC(P,B,Lmax)
% [E,F]=PBC2EFC(P,B,Lmax)
%
% Transforms the vector spherical harmonic coefficients from the PBC basis
% (classical) to the EFC basis (up/downward continuable). 
%
% Two different modes: 
% EITHER the P and B coefficients are given as Blmcosi, Blmcosi;
% then the returning E F will also be in the lmcosi format (Elmcosi, 
% Flmcosi),
% OR P and B are given in either the ADDMOUT or ADDMON format 
% (both work because the factor depends on the degree, not on the order).
% In this case the returning coefficients will be given in the same format
%
% Remarks: 1) C is not required, because it is not transformed
%          2) E has degrees 0,1,... F has degrees 1,... 
%
% INPUT: 
%
% P, B  Radial and Gradient projection component coefficients. Either in
%       the lmcosi or in the ADDMOUT or ADDMON format. If given in the 
%       ADDMOUT or ADDMON format, Lmax is required and all coefficients are
%       needed.
% Lmax  Maximum degree. Only needed if the coefficients are given in either
%       the ADDMON or ADDMOUT format.
%
% OUTPUT:
%
% E, F  coefficients for the up/downward continuable vector spherical
%       harmonics (eigenfunctions of the vectorial Beltrami operator). 
%       Either given as lmcosi or vector of coefficients (depending on the 
%       input).
%
% EXAMPLE:
%
% PBC2EFC('demo1')  Generate random P B. Transform to E F. Check if lmcosi
% and vector produce the same. Look at P B and E F norm to show that the
% transformation is orthogonal.
%
% COMMENT:
%
% This is done using a for loop. I guess it would be faster if the loop is
% replaced by a vector multiplication.
%
% Last modified by plattner-at-alumni.ethz.ch, 08/20/2012
%
% See also EFC2PBC

defval('Lmax',-1);

if ~ischar(P)

% Decision if P and B are given as lmcosi
if (size(P,2)==4)&&(size(B,2)==4)&&(Lmax==-1)
    % The lmcosi case 
    
    % Remove the P_00 coefficient from the list to make things easier
    Poo=P(1,:);
    P=P(2:end,:);
    
    E=nan(size(P)); F=nan(size(B));
    % Copy the l and m values
    E(:,1)=P(:,1); E(:,2)=P(:,2);
    F(:,1)=B(:,1); F(:,2)=B(:,2);
    
    for i=1:size(B,1)
        L=B(i,1);
        Efac1= sqrt((L+1)/(2*L+1));
        Efac2=-sqrt((L  )/(2*L+1));
        Ffac1=-Efac2;
        Ffac2= Efac1;

        E(i,3)=Efac1*P(i,3) + Efac2*B(i,3); 
        E(i,4)=Efac1*P(i,4) + Efac2*B(i,4);

        F(i,3)=Ffac1*P(i,3) + Ffac2*B(i,3); 
        F(i,4)=Ffac1*P(i,4) + Ffac2*B(i,4); 

    end
    
    % And now put the P_00 coefficient into E
    E=[Poo;E];
    
elseif (size(P,1)==(Lmax+1)^2) && (size(B,1)==(Lmax+1)^2-1)
    % The list case    
    
    % Remove the P_00 coefficient to make things easier
    Poo=P(1,:);
    P=P(2:end,:);
    
    [~,~,~,~,~,~,~,bigl]=addmon(Lmax);
    bigl=bigl(2:end);
    E=nan(size(P)); F=nan(size(B));
    for i=1:(Lmax+1)^2-1
        L=bigl(i);
        Efac1= sqrt((L+1)/(2*L+1));
        Efac2=-sqrt((L  )/(2*L+1));
        Ffac1=-Efac2;
        Ffac2= Efac1;

        E(i,:)=Efac1*P(i,:) + Efac2*B(i,:); 
        F(i,:)=Ffac1*P(i,:) + Ffac2*B(i,:);
    end
    
    % And now put the P_00 coefficient into E
    E=[Poo;E];
    
else
    error('Something is wrong with the input. Or you use vector and forgot to provide Lmax')
end


elseif strcmp(P,'demo1')
   % Create P and B randomly and then transform it using both, lmcosi and 
   % vector formulation
   Lmax=300;
   Plmcosi=plm2rnd(Lmax);
  
   Blmcosi=plm2rnd(Lmax);
   Blmcosi=Blmcosi(2:end,:);
   
   % First use lmcosi
   [Elmcosi,Flmcosi]=PBC2EFC(Plmcosi,Blmcosi);
   
   vec=blmclm2coef(Elmcosi(2:end,:),Flmcosi);
   vec=vec(:);
   E1=[Elmcosi(1,3);vec(1:end/2)]; F1=vec(end/2+1:end);
   
   % Then use vector formulation
   PBvec=blmclm2coef(Plmcosi(2:end,:),Blmcosi);
   PBvec=PBvec(:);
   P=[Plmcosi(1,3);PBvec(1:end/2)];
   B=PBvec(end/2+1:end);
   
   [E2,F2]=PBC2EFC(P,B,Lmax);
   
   EF1=[E1;F1];
   EF2=[E2;F2];
   PB=[P;B];
   fprintf('Norm difference between the methods is %g\n',norm(EF1-EF2));
   fprintf('Norm of PB coef is %g, norm of EF coef is %g\n',...
       norm(EF1),norm(PB));
 
   
   
   
end
    
function [V,C,Vm,Cm]=vectansdwcap(TH,L,m)
% [V,C,Vm,Cm]=vectansdwcap(TH,L,m)
%
% Spherical-harmonic localization to a single fixed-order spherical polar
% cap, bandlimited and optimally spatially concentrated solutions. For the
% tangential component of the vector Slepian functions.
%
% INPUT:
%
% TH    angular extent of the spherical cap (semi-opening angle) in degrees
%       or two angles for ring between
% L     maximum spherical-harmonic degree
% m     spherical-harmonic order
%
% OUTPUT:
%
% V     eigenvalues for positive (or zero) m
% C     eigenvectors = coefficients for the tangential vector spherical
%       harmonics Blm and Clm for positive (or zero) m
% Vm    eigenvalues for negative m
% Cm    eigenvectors = coefficients for the tangential vector spherical
%       harmonics Blm and Clm for negative m
%
% Last modified by plattner-at-alumni.ethz.ch, 5/8/2017


dirname=fullfile(getenv('IFILES'),'VECTANSDWCAP');

fnpl=fullfile(dirname,sprintf(...
      'VECTANSDW-%f-%f-%i-%i-%i.mat',min(TH),max(TH),min(L),max(L),m));

if exist(fnpl,'file')==2 
    load(fnpl)
    disp(sprintf('%s loaded by CAPVECTORSLEPIAN',fnpl))
else 
    if length(TH)==2
        [Mm1,Mmm1]=kerneltancapm(max(TH),L,m);
        [Mm2,Mmm2]=kerneltancapm(min(TH),L,m);
        Mm=Mm1-Mm2;
        Mmm=Mmm1-Mmm2;
    else
        [Mm,Mmm]=kerneltancapm(TH,L,m);
    end
    % Positive m
    [C,V]=eig(Mm);
    [V,isrt]=sort(sum(V,1),'descend');
    C=C(:,isrt);
    % Negative m
    if m~=0
        [Cm,Vm]=eig(Mmm);
        [Vm,isrtm]=sort(sum(Vm,1),'descend');
        Cm=Cm(:,isrtm);
    else
        Cm=C; Vm=V;
    end
    save(fnpl,'C','V','Cm','Vm');
end


  
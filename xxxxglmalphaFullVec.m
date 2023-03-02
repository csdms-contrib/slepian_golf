function [G,V] = glmalphaFullVec(Lmax,dom,PorE,rotcoord,J,OutOrOn)
% [G,V] = glmalphaFullVec(Lmax,dom,PorE,J,rotcoord,OutOrOn)
%
% Coefficients for full vector Slepian functions in either the Plm Blm Clm,
% or Elm Flm Clm basis.
%
% INPUT:
% Lmax          Maximum spherical harmonic degree
% dom           Region name or semi-opening angle in degrees
% PorE          Which basis? EFC=1 or true, PBC=0 or false
% rotcoord     If rotated spherical cap: [lon, colat] of cap center, both in degrees.
%                   otherwise (if named region or polar cap): set it to []
% J                 If rotated spherical cap: How many Slepian functions to
%                   write out? 
% OutOrOn    ADDMOUT (1 or true) format or ADDMON (0 or false) format?
% 
% OUTPUT:
% G          Spherical harmonic coefficients to build Slepian functions
%            in ADDMOUT order (unless asking for ADDMON)
% V          Concentration values for the Slepian functions
% 
% Last modified by plattner-at-alumni.ethz.ch, 2023/3/2

defval('PorE', false)
defval('rotcoord', [])
defval('J', 3*(Lmax+1)^2-2)
defval('OutOrOn', false)

% Create radial kernel and solve for Slepian functions
[Gr,Vr] = glmalpha(dom,Lmax);

% Create tangential kernel and solve for Slepian functions
[Gt,Vt] = vectanglmalpha(dom,Lmax);

% Combine and order the Slepian functions
G=zeros(length(Gr)+length(Gt));
G(1:length(Gr),1:length(Gr))=Gr;
G(length(Gr)+1:end,length(Gr)+1:end)=Gt;
V=[Vr(:);Vt(:)];

[V,isrt] = sort(V,'descend');
G = G(:,isrt);


% If rotated spherical cap, rotate the best J Slepian functions
if ~isempty(rotcoord)
    fname=fullfile(getenv('IFILES'),'FULLVECROTJ',...
        sprintf('glmalphaFullVecRotJ-%g-%i-%g-%g-%i.mat',...
        dom,Lmax,rotcoord(1),rotcoord(2),J));

    if exist(fname,'file')==2
        load(fname)
        disp(sprintf('Loading %s',fname))
    else
        % Calculate
        Grot = nan(size(G,1),J);
        %parfor  j=1:J
        for  j=1:J
            % Rotating each P,B,C component individually. This should work
            % because of the linearity of the rotation.
            prng = 1:length(Gr);
            brng = length(Gr)+1: length(Gr)+length(Gt)/2;
            crng = length(Gr)+length(Gt)/2+1: length(Gr)+length(Gt);
            lmcosiP = coef2lmcosi( G(prng,j), 1);
            lmcosiB = coef2lmcosi(  [0; G(brng,j)] ,1);
            lmcosiC = coef2lmcosi(  [0; G(crng,j)] ,1);

            lmcosiProt=plm2rot(lmcosiP, 0, -rotcoord(2),-rotcoord(1));
            lmcosiBrot=plm2rot(lmcosiB, 0, -rotcoord(2),-rotcoord(1));
            lmcosiCrot=plm2rot(lmcosiC, 0, -rotcoord(2),-rotcoord(1));

            Prot = lmcosi2coef(lmcosiProt, 1);
            Brot = lmcosi2coef(lmcosiBrot, 1);
            Brot = Brot(2:end);
            Crot = lmcosi2coef(lmcosiCrot, 1);
            Crot = Crot(2:end);
            Grot(:,j) = [Prot;Brot;Crot];

        end

    G = Grot;
    V = V(1:J);
    % Now save
    try
        % If you are running Matlab
        save(fname,'G','V','-v7.3')
    catch
        % If you are running octave
        save(fname,'G','V')
    end

    end % Done with calculation

end

% If EFC basis, transform the J best Slepian functions.
if PorE
    [E,F] = PBC2EFC( G(1:(Lmax+1)^2,:), G((Lmax+1)^2+1:2*(Lmax+1)^2-1,:), Lmax);

    G(1:(Lmax+1)^2,:) = E;
    G((Lmax+1)^2+1:2*(Lmax+1)^2-1,:) = F;
end


% OutOrOn
if ~OutOrOn
    G = out2on(G,Lmax);
end











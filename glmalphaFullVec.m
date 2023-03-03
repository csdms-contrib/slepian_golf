function [G,V] = glmalphaFullVec(Lmax,dom,PorE,rotcoord,J,OnOrOut)
  % [G,V] = glmalphaFullVec(Lmax,dom,PorE,J,rotcoord,OnOrOut)
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
  % J                 If rotated spherical cap or large Lmax: How many Slepian functions to
  %                   write out? 
  % OnOrOout    ADDMOUT (1 or true) format or ADDMON (0 or false) format?
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
  defval('OnOrOut', false)

  
  % Create radial kernel and solve for Slepian functions
  [Gr,Vr] = glmalpha(dom,Lmax);
Gr = sparse(Gr);

  % Create tangential kernel and solve for Slepian functions
  [Gt,Vt] = vectanglmalpha(dom,Lmax);
Gt = sparse(Gt);


  % Combine and order the Slepian functions
  G=sparse(length(Gr)+length(Gt), length(Gr)+length(Gt));
  G(1:length(Gr),1:length(Gr))=Gr;
  G(length(Gr)+1:end,length(Gr)+1:end)=Gt;
  V=[Vr(:);Vt(:)];
  clear Gr; clear Gt; clear Vr; clear Vt; 


  [V,isrt] = sort(V,'descend');
  G = G(:,isrt);
  % Now keep only the J best
  G = G(:,1:J);

  if isempty(rotcoord) | isstr(dom)
      % OnOrOut
      if ~OnOrOut
          G(1:(Lmax+1)^2, :) = out2on( G(1:(Lmax+1)^2,:), Lmax);
          tmp = out2on( [zeros(1,size(G,2));   G( (Lmax+1)^2+1:2*(Lmax+1)^2-1, :) ], Lmax);
          G( (Lmax+1)^2+1:2*(Lmax+1)^2-1, :) = tmp(2:end,:);
          tmp = out2on( [zeros(1,size(G,2));  G(2*(Lmax+1)^2:3*(Lmax+1)^2-2, :) ], Lmax);
          G( 2*(Lmax+1)^2:3*(Lmax+1)^2-2, :) = tmp(2:end,:);
          clear tmp
      end

  % If rotated spherical cap, rotate the best J Slepian functions
  else
    fname=fullfile(getenv('IFILES'),'FULLVECROTJ',...
                   sprintf('glmalphaFullVecRotJ-%g-%i-%g-%g-%i.mat',...
                           dom,Lmax,rotcoord(1),rotcoord(2),J));

    if exist(fname,'file')==2
      load(fname)
      disp(sprintf('Loading %s',fname))
      % This is stored in ADDMON
    else
      % Calculate

      %% Need to free up memory... Do exactly as in gradvecglmalphauptoJp.m
      G = sparse(G(:,1:J));

      % out2on for G
      G(1:(Lmax+1)^2, :) = out2on( G(1:(Lmax+1)^2,:), Lmax);
      tmp = out2on( [zeros(1,size(G,2));   G( (Lmax+1)^2+1:2*(Lmax+1)^2-1, :) ], Lmax);
      G( (Lmax+1)^2+1:2*(Lmax+1)^2-1, :) = tmp(2:end,:);
      tmp = out2on( [zeros(1,size(G,2));  G(2*(Lmax+1)^2:3*(Lmax+1)^2-2, :) ], Lmax);
      G( 2*(Lmax+1)^2:3*(Lmax+1)^2-2, :) = tmp(2:end,:);
      clear tmp

      V = V(1:J);
      Grot =sparse(size(G,1),size(G,2));
      [demsz,delsz,mz,lmc,mzin,mzo]=addmon(Lmax);
      prng = 1:(Lmax+1)^2;
      brng = (Lmax+1)^2+1:2*(Lmax+1)^2-1;
      crng = 2*(Lmax+1)^2:3*(Lmax+1)^2-2;
      
      %parfor  j=1:J
      for  j=1:J
        % Rotating each P,B,C component individually.    
        lmcosiP = full([delsz,demsz,reshape(insert(G(prng,j),0,mzin),2,length(demsz))']);
        lmcosiProt=plm2rot(lmcosiP, 0, -rotcoord(2),-rotcoord(1));
        %clear lmcosiP;
        cosivec=reshape(lmcosiProt(:,3:4),1,2*length(demsz));
        Prot=cosivec(mzo);
        %clear lmcosiProt;
        
        lmcosiB = full([delsz,demsz,reshape(insert([0;G(brng,j)],0,mzin),2,length(demsz))']);
        lmcosiBrot=plm2rot(lmcosiB, 0, -rotcoord(2),-rotcoord(1));
        %clear lmcosiB;
        cosivec=reshape(lmcosiBrot(:,3:4),1,2*length(demsz));
        Brot=cosivec(mzo);
        %clear lmcosiBrot;
        Brot = Brot(2:end);
        
        lmcosiC = full([delsz,demsz,reshape(insert([0;G(crng,j)],0,mzin),2,length(demsz))']);
        lmcosiCrot=plm2rot(lmcosiC, 0, -rotcoord(2),-rotcoord(1));
        %clear lmcosiC;   
        cosivec=reshape(lmcosiCrot(:,3:4),1,2*length(demsz));
        Crot=cosivec(mzo);
        %clear lmcosiCrot;
        Crot = Crot(2:end);
        
        Grot(:,j) = [Prot(:);Brot(:);Crot(:)];
        %clear Prot; clear Brot, clear Crot;        
      end

      G = Grot;
      clear Grot;
      % Now save
      try
        % If you are running Matlab
        save(fname,'G','V','-v7.3')
      catch
        % If you are running octave
        save(fname,'G','V')
      end

    end % Done with calculation


    %%%%%% need to do onorout. G is calculated / stored in addmon, because of rotation.
    %%%%%% If the user wants addmout, need to take care of this.
    if OnOrOut
        [~,~,~,~,~,~,~,~,rinm]=addmon(Lmax);
        G(1:(Lmax+1)^2,:) = G(rinm,:);
        G((Lmax+1)^2+1:2*(Lmax+1)^2-1, :) = G((Lmax+1)^2+rinm(2:end)-1,:);
        G(2*(Lmax+1)^2:3*(Lmax+1)^2-2, :) = G(2*(Lmax+1)^2-1+rinm(2:end)-1,:);
    end

  end

  % If EFC basis, transform the J best Slepian functions.
  if PorE
    [E,F] = PBC2EFC( G(1:(Lmax+1)^2,:), G((Lmax+1)^2+1:2*(Lmax+1)^2-1,:), Lmax);

    G(1:(Lmax+1)^2,:) = E;
    G((Lmax+1)^2+1:2*(Lmax+1)^2-1,:) = F;
  end














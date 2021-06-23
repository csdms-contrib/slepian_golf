function varargout=inoutgradvecglmalpha(TH,Lin,Lout,srt)
% [G,V]=inoutgradvecglmalpha(TH,Lin,Lout,srt)
%
% Construction of mixed internal- and external source gradient vector Slepian 
% functions from the the Elm and Flm with two different max
% spherical-harmonic degrees
%
% INPUT:
%
% TH    Named region or polar cap opening angle
% Lin   maximum spherical harmonic degree for inner sources
% Lout  maximum spherical harmonic degree for outer sources
% srt       Should the Slepian functions be sorted? [default = 1 = yes]
%
% OUTPUT:
%
% G     Matrix containing the slepian coefficients for linear combinations
%       of the Elm and Flm. The first (Lin+1)^2 coefficients are for the
%       Elm, the last (Lout+1)^2-1 coefficients are for the Flm
% V     suitability (eigen) values
%
% Last modified by plattner-at-alumni.ethz.ch, 05/05/2017

defval('anti',0)
defval('srt',1)

% First check if already calculated. Make the name:
if ~ischar(TH)     
    if length(TH)==1 % POLAR CAPS  
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHA',...
		  sprintf('inoutgradvecglmalpha-%g-%i-%i-%i.mat',TH,max(Lin),min(Lin),Lout));
    elseif length(TH)==2          
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHA',...
		  sprintf('inoutgradvecglmalpha-%g-%g-%i-%i-%i.mat',...
          max(TH),min(TH),max(Lin),min(Lin),Lout));   
    else
        error('Provide one or two spherical cap opening angles')
    end
else % GEOGRAPHICAL REGIONS and XY REGIONS
  if length(Lin)>1
    warning('Bandpass for internal field not yet implemented. Working with max(Lin)');
    Lin=max(Lin);
  end
 % This is in case we give the region as a list of lon lat coordinates
  if ischar(TH) % Geographic (keep the string)
    h=TH;
  else % Coordinates (make a hash)
    h=hash(TH,'sha1');
  end   
  fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHA',...
		 sprintf('inoutgradvecglmalpha-%s-%i-%i.mat',h,Lin,Lout));         
end


if exist(fname,'file')==2
    load(fname)
    fprintf('Loading %s\n',fname)
else       
    % If possible, do parallel
    try 
        parpool
    %catch
    %    disp('Parpool already open')
    end

    % Find row indices into G belonging to the orders
    [EM,EL,mz,blkm]=addmout(max(Lin)); 

    % For GEOGRAPHICAL REGIONS or XY REGIONS
    if ischar(TH) || length(TH)>2    
        
        % Calculate the ElmElm localization kernel part
        try
            KEE=kernelep(Lin,TH);
        catch
            KEE=kernele(Lin,TH); 
        end
    
        % Calculate the FlmFlm localization kernel part
        try
            KFF=kernelfp(Lout,TH);
        catch
            KFF=kernelf(Lout,TH); 
        end
  
        % Now the mixed component. This should be quick, only loading
        KEF=kernelmixefp(Lin,Lout,TH);
      
        % Now put them all together
        K=[KEE KEF;KEF' KFF];
        % And avoid asymmetry caused by numerical noise
        K=(K+K')/2;
        
        % Calculate the eigenvectors / eigenvalues
        [G,V]=eig(K);

        % Sort eigenvectors for descending eigenvalues
        [V,isrt]=sort(sum(real(V),1));
        V=fliplr(V);
        G=G(:,fliplr(isrt));  
      
        % Reorder eigenvector coefficient ordering from ADDMON to
        % ADDMOUT
        [~,~,~,~,~,~,~,~,Rin,~]=addmon(Lin);
        [~,~,~,~,~,~,~,~,Rout,~]=addmon(Lout);
        
        resort=[Rin;Rout];
        % To do addmon to addmout: include the Lout=0 row
        Gresort=[G(1:(Lin+1)^2,:);nan(1,size(G,2));G((Lin+1)^2+1:end,:)];
        % Then resort the individual parts
        Gresort=Gresort(resort,:);
        % Then remove the Lout=0 row
        G=[Gresort(1:(Lin+1)^2,:);Gresort((Lin+1)^2+2:end,:)];
        
        N=sum(V);
        EL=[];
        EM=[];
        
        save(fname,'G','V','EL','EM','N')  
        
    else

      % Polar caps
      % Plattner 5/5/2017:
      % For now I'm using the strategy: Set up the entire max(Lin) part but
      % the first min(Lin)^2+1 entries will be zero.    

      maxLin=max(Lin);
      
      deME=addmout(maxLin);
      deMF=addmout(Lout);

      % In the alpha vector we need to include the correct sizes of the
      % matrices
      mvec=0:max(maxLin,Lout);
      if length(Lin)==2
        sizE=maxLin+1-max(mvec,min(Lin));
      else
        sizE=subplus(maxLin+1-mvec);
      end
      sizF=subplus(Lout+1-max(mvec,1));    

      % The vectors sizE and sizF contain the sizes of the matrices for the E
      % components and the F components. Of course we will have them both.
      siztot=sizE+sizF;

      % The entries in alpha show the beginning (alpha(i)) and end 
      % (alpha(i+1))-1 of each block i. The first block, for m=0, only occurs
      % once. then for m=-1 m=1, etc we always have the positive and negative
      % m. Therefore we start with index 1, then the size of the m=0 block, 
      % then all other sizes twice. That's what gamini does. Then of course 
      % we need the cumulative index to know where in the big matrix we want
      % to put things.
    
      if length(Lin==2)
        % We need to + and - m, but m=0 is only once
        sizindices=[1 gamini(2:max(maxLin,Lout)+1,2)];
        % This could be done in a direct way without the for loop, but it's not
        % expensive to do this here..
        alpha=[1;nan(length(sizindices)-1,1)];
        for i=2:length(alpha)
            alpha(i)=alpha(i-1)+siztot(sizindices(i-1));
        end
        % The last one has the same number of Ls as the 
        % previous to last one because of m, -m 
        alpha=[alpha;alpha(end)+(alpha(end)-alpha(end-1))];
    else
        alpha=cumsum([1 siztot(1) gamini(siztot(2:end),2) ]);
    end    

    % The redistribution is not even as complicated as it seems. The alpha
    % intervalls need to be exactly the size of the number of eigenvectors.
    % The right location with the ms is taken care of by EM==m 
    % (might need to make this faster?)
    
    % To simplify things and avoid mistakes: Treat the Flm coefficients the
    % same way as the Elm coefficients and just remove the L=0 part if
    % necessary.
    
    % Initialize matrices  
    if length(Lin)==2
        GE=sparse((maxLin+1)^2,(maxLin+1)^2 - (min(Lin))^2 + (Lout+1)^2-1);
        % Treat Flm like Elm and remove the zero part later
        GF=sparse((Lout+1)^2,  (maxLin+1)^2 - (min(Lin))^2 + (Lout+1)^2-1);     
    else
        GE=sparse((maxLin+1)^2,(Lin+1)^2 + (Lout+1)^2-1);
        % Treat Flm like Elm and remove the zero part later
        GF=sparse((Lout+1)^2,  (Lin+1)^2 + (Lout+1)^2-1); 
    end
    V=zeros(1,(maxLin+1)^2 - (min(Lin)+1)^2 + (Lout+1)^2-1);
    disp('Calculating in parallel mode')

    % Now the same as in inoutgradvecglmalphaup: Calculate the individual
    % solutions for the m and put them back in the right place
    parfor mm=1:max(maxLin,Lout)+1       
    %for mm=1:max(maxLin,Lout)+1   
            m=mm-1;            
            [Vpp,Cp]=inoutgradvecsdwcap(TH,Lin,Lout,m);
            Vp{mm}=Vpp;
            sizE=subplus(maxLin+1-m);
            %sizF=subplus(Lout+1-max(m,1));
            CE{mm}=Cp(1:sizE,:);            
            CF{mm}=Cp(sizE+1:end,:);
    end 

    % To make distribution a bit simpler, add the L=0 row to CF.
    % This is only necessary for m=0.
    CF{1}=[zeros(1,size(CF{1},2)); CF{1}];
    
    % Distribute this at the right point in the huge matrix
    for m=0:max(maxLin,Lout)   
        if m>0
            % Here you supply the negative orders
            GE(deME==-m,alpha(2*m):alpha(2*m+1)-1)=CE{m+1};  
            GF(deMF==-m,alpha(2*m):alpha(2*m+1)-1)=CF{m+1};  
            V(alpha(2*m):alpha(2*m+1)-1)=Vp{m+1};
        end
        % Duplicate for the positive order in case the region is axisymmetric
            GE(deME==m,alpha(2*m+1):alpha(2*m+2)-1)=CE{m+1};     
            GF(deMF==m,alpha(2*m+1):alpha(2*m+2)-1)=CF{m+1};  
            V(alpha(2*m+1):alpha(2*m+2)-1)=Vp{m+1};   
    end           

    GF=GF(2:end,:);
    G=[GE;GF];

    if srt
        [V,isrt]=sort(V,'descend');
        % Now remove the L=0 part from GF. It should be a zero row anyhow
        % Uncomment this if you want to test it
        % fprintf('Norm of L=0 in GF is %g\n',norm(GF(1,:)))        
        G=G(:,isrt);
    end
        
    try
    	% Matlab
    	save(fname,'G','V','-v7.3')
    catch
      % Octave
    	save(fname,'G','V')
    end
    
      
    end % end of polar caps
    
end % end of calculation if not yet available

% Provide output
varns={G,V};
varargout=varns(1:nargout);   
         

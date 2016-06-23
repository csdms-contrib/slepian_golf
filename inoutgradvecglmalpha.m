function varargout=inoutgradvecglmalpha(TH,Lin,Lout)
% [G,V]=inoutgradvecglmalpha(TH,Lin,Lout)
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
%
% OUTPUT:
%
% G     Matrix containing the slepian coefficients for linear combinations
%       of the Elm and Flm. The first (Lin+1)^2 coefficients are for the
%       Elm, the last (Lout+1)^2-1 coefficients are for the Flm
% V     suitability (eigen) values
%
% Last modified by plattner-at-alumni.ethz.ch, 04/15/2015



% First check if already calculated. Make the name:
if ~ischar(TH)     
    if length(TH)==1 % POLAR CAPS  
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHA',...
		  sprintf('inoutgradvecglmalpha-%g-%i-%i.mat',TH,Lin,Lout));
    elseif length(TH)==2          
        fname=fullfile(getenv('IFILES'),'INOUTGRADVECGLMALPHA',...
		  sprintf('inoutgradvecglmalpha-%g-%g-%i-%i.mat',...
          max(TH),min(TH),Lin,Lout));   
    else
        error('Provide one or two spherical cap opening angles')
    end
else % GEOGRAPHICAL REGIONS and XY REGIONS
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
        matlabpool open
    %catch
    %    disp('Matlabpool already open')
    end
    
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
        error('Not implemented for spherical caps')
    end
      
end

% Provide output
varns={G,V};
varargout=varns(1:nargout);   
         

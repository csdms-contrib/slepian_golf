function varargout=outgradvecglmalpha(TH,L)
% [G,V]=outgradvecglmalpha(TH,L)
%
% Construction of external-field gradient vector Slepian functions from the
% Flm.
%
% INPUT:
%
% TH    Named region or polar cap opening angle
% L     maximum spherical harmonic degree
%
% OUTPUT:
%
% G     Matrix containing the slepian coefficients for linear combinations
%       of the Flm
% V     concentration values
%
% Last modified by plattner-at-alumni.ethz.ch, 02/24/2015

% Figure out if it's lowpass or bandpass
lp=length(L)==1;
bp=length(L)==2;
maxL=max(L);

% The spherical harmonic dimension
ldim=(L(2-lp)+1)^2-bp*L(1)^2;
% if we start with L=0 we need to remove L=0
if lp
    ldim=(L+1)^2-1;
end

% First check if already calculated. Make the name:
if ~isstr(TH) && length(TH)==1 % POLAR CAPS  
    if lp
      fname=fullfile(getenv('IFILES'),'OUTGRADVECGLMALPHA',...
		     sprintf('outgradvecglmalpha-%g-%i.mat',TH,L));
    elseif bp
      fname=fullfile(getenv('IFILES'),'OUTGRADVECGLMALPHA',...
		     sprintf('outgradvecglmalphabl-%g-%i-%i.mat',TH,L(1),L(2)));
    else
      error('The degree range is either one or two numbers')       
    end
    % Initialize ordering matrices
    %MTAP=repmat(0,1,ldim);
    %IMTAP=repmat(0,1,ldim);
    
else % GEOGRAPHICAL REGIONS and XY REGIONS
    % This is in case we give the region as a list of lon lat coordinates
    if isstr(TH) % Geographic (keep the string)
      h=TH;
    else % Coordinates (make a hash)
      h=hash(TH,'sha1');
    end
    if lp
      fname=fullfile(getenv('IFILES'),'OUTGRADVECGLMALPHA',...
		     sprintf('outgradvecglmalpha-%s-%i.mat',h,L));
    elseif bp
      fname=fullfile(getenv('IFILES'),'OUTGRADVECGLMALPHA',...
		     sprintf('outgradvecglmalphabl-%s-%i-%i.mat',h,L(1),L(2)));
    else
     error('The degree range is either one or two numbers')       
    end  
    
end
  
  
  
% We have the name. Now see if it is already calculated  
if exist(fname,'file')==2
    load(fname)
    disp(sprintf('Loading %s',fname))
else  
    % Do all the calcualtions!!
    % Find row indices into G belonging to the orders
    [EM,EL,mz,blkm]=addmout(maxL);
    % Find increasing column index; that's how many belong to this order
    % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
    % The middle bit is twice for every nonzero order missing
    % alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
    %   		gamini(L(2-lp)-bp*(L(1)-1),bp*2*(L(1)-1)) ...
    %   		gamini(L(2-lp)-bp*(L(1)-1):-1:1,2)]);
    % This should be the same for L and [0 L]
    alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
        gamini(L(2-lp)-bp*(L(1)-1),bp*2*L(1)) ...
        gamini(L(2-lp)-bp*L(1):-1:1,2)]);
    % To go from alpha for L=0 to alpha starting with L=1, we need to make
    % the first block by one smaller and keep the following ones the same
    % size but shifted by one location.
    alpha(1)=alpha(1)+1;
    alpha=alpha-1;
    % And also remove the L=0 entry from the list of m values. This is
    % easy:
    EM=EM(2:end);
    
    % For GEOGRAPHICAL REGIONS or XY REGIONS
    if ischar(TH) || length(TH)>1
        if bp
            error('Bandpass geographical tapers are not ready yet')
        end
        % Calculates the localization kernel for this domain
        Klmlmp=kernelfp(L,TH);

        % Calculate the eigenvectors / eigenvalues
        [G,V]=eig(Klmlmp);

        % Sort eigenvectors for descending eigenvalues
        [V,isrt]=sort(sum(real(V),1));
        V=fliplr(V);
        G=G(:,fliplr(isrt));  

        % Reorder eigenvector coefficient ordering from ADDMON to
        % ADDMOUT
        [a,b,c,d,e,f,ems,els,R1,R2]=addmon(L);
        % Remember, we don't have degree 0. But degree 0 does not
        % switch its location between ADDMON and ADDMOUT. It's always
        % the first
        % Must make G one larger
        G=[zeros(1,size(G,2));G];
        G=G(R1,:);
        G=G(2:end,:);

        % Calculate Shannon number
        N=sum(V);

        save(fname,'G','V','EL','EM','N')  

    else
        % For AXISYMMETRIC REGIONS
        % Initialize matrices
        G=sparse((maxL+1)^2-1,ldim);%repmat(0,(maxL+1)^2,ldim);
        V=zeros(1,ldim); 

        for m=0:maxL
            [~,Vp,~,~,C]=outgradvecsdwcap(TH,L,m,0,-1);
            % Distribute this at the right point in the huge matrix
            if m>0
                % Here you supply the negative orders
                % Maybe can make this faster sometime by building
                % individual small matrices and then directly assemble
                % them in a sparse way.
                G(EM==-m,alpha(2*m):alpha(2*m+1)-1)=C;
                V(alpha(2*m):alpha(2*m+1)-1)=Vp;
                %MTAP(alpha(2*m):alpha(2*m+1)-1)=-m;
                % It's all neatly ordered here, 
                % downgoing within every order
                %IMTAP(alpha(2*m):alpha(2*m+1)-1)=1:length(Vp);                    
            end                
            % Duplicate for the positive order in case the region is axisymmetric
            % Maybe can make this faster sometime by building
            % individual small matrices and then directly assemble
            % them in a sparse way.
            G(EM==m,alpha(2*m+1):alpha(2*m+2)-1)=C;
            V(alpha(2*m+1):alpha(2*m+2)-1)=Vp;
            %MTAP(alpha(2*m+1):alpha(2*m+2)-1)=m;
            % It's all neatly ordered here, downgoing within every order
            %IMTAP(alpha(2*m+1):alpha(2*m+2)-1)=1:length(Vp);
        end
        
        if srt
            [V,isrt]=sort(V,'descend');
            G=G(:,isrt);
        end
        
        save(fname,'G','V','EL','EM')

    end
end
    
    
% Provide output
varns={G,V};
varargout=varns(1:nargout);                
                    

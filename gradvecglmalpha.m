function varargout=gradvecglmalpha(TH,L,srt)
% [G,V]=gradvecglmalpha(TH,Lmax,srt)
% 
% Construction of internal-field gradient vector Slepian functions
%
% INPUT:
%
% TH    Named region or cpolar cap opening angle
% L     maximum spherical harmonic degree
% srt   sorted output?
%
% OUTPUT:
%
% G     Matrix containing the slepian coefficients for linear combinations
%       of the Elm
% V     concentration values
%
% Last modified by plattner-at-alumni.ethz.ch, 02/26/2015

defval('srt',1)
defval('anti',0)

% Figure out if it's lowpass or bandpass
lp=length(L)==1;
bp=length(L)==2;
maxL=max(L);

% The spherical harmonic dimension
ldim=(L(2-lp)+1)^2-bp*L(1)^2;

% First check if already calculated
  if ~isstr(TH) && length(TH)==1 % POLAR CAPS
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHA',...
		     sprintf('gradvecglmalpha-%g-%i.mat',TH,L));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHA',...
		     sprintf('gradvecglmalphabl-%g-%i-%i.mat',TH,L(1),L(2)));
    else
      error('The degree range is either one or two numbers')       
    end
    
       % Initialize ordering matrices
    %MTAP=repmat(0,1,ldim);
    %IMTAP=repmat(0,1,ldim);
  else % GEOGRAPHICAL REGIONS and XY REGIONS
    defval('sord',10) % SPLINING SMOOTHNESS
    % We'll put in a Shannon number based on the area only, not based on
    % an actual sum of the eigenvalues
    defval('J',ldim)
    % Note the next line, though we can change our minds
    defval('J',ldim*spharea(TH))
    if isstr(TH) % Geographic (keep the string)
      h=TH;
    else % Coordinates (make a hash)
      h=hash(TH,'sha1');
    end
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHA',...
		     sprintf('gradvecglmalpha-%s-%i-%i.mat',h,L,J));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHA',...
		     sprintf('gradvecglmalphabl-%s-%i-%i-%i.mat',h,L(1),L(2),J));
    else
     error('The degree range is either one or two numbers')       
    end
    defval('GM2AL',NaN) % If not, calculate order per taper
    defval('MTAP',NaN) % If not, calculate order per taper
    defval('IMTAP',NaN) % And rank ordering within that taper
    defval('xver',0) % For excessive verification of the geographical case
  end
  
if exist(fname,'file')==2
  load(fname)
  disp(sprintf('Loading %s',fname))
else
  % Initialize matrices
  % plattner-at-alumni.ethz.ch, 2/11/2015:
  % This is now moved to the individual cases because generic regions will
  % reqire full matrices and polar caps require sparse matrices  
  %G=repmat(0,(maxL+1)^2,ldim);
  %V=repmat(0,1,ldim);
  
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
      
  % For GEOGRAPHICAL REGIONS or XY REGIONS
  if ischar(TH) || length(TH)>1
    % Initialize matrices
    %G=zeros(0,(maxL+1)^2,ldim);
    %V=zeros(0,1,ldim);   
    if bp
      error('Bandpass geographical tapers are not ready yet')
    end
    % Calculates the localization kernel for this domain
      Klmlmp=kernelep(L,TH);

    
    if anti==1
      % Get the complimentary region
      Klmlmp=eye(size(Klmlmp))-Klmlmp;
    end
    
    
    % Calculates the eigenfunctions/values for this localization problem
    [G,V]=eig(Klmlmp);
    [V,isrt]=sort(sum(real(V),1));
    V=fliplr(V);
    G=G(:,fliplr(isrt));
    
    [a,b,c,d,e,f,ems,els,R1,R2]=addmon(L);
    % This indexes the orders of G back as 0 -101 -2-1012 etc
    G=G(R1,:);
    % Check indexing
%     difer(els(R1)-EL,[],[],mesg)
%     difer(ems(R1)-EM,[],[],mesg)
    
    % Calculate Shannon number and compare this with the theory
    N=sum(V);
    G=G(:,1:J);
    V=V(1:J);
    save(fname,'G','V','EL','EM','N')  
    
  else
        % For AXISYMMETRIC REGIONS
        % Initialize matrices
        G=sparse((maxL+1)^2,ldim);%repmat(0,(maxL+1)^2,ldim);
        V=zeros(1,ldim);
        for m=0:maxL
            [E,Vp,Np,th,C]=gradvecsdwcap(TH,L,m,0,-1);
             % Distribute this at the right point in the huge matrix
        if m>0
            % Here you supply the negative orders
            G(EM==-m,alpha(2*m):alpha(2*m+1)-1)=C;
            V(alpha(2*m):alpha(2*m+1)-1)=Vp;
            %MTAP(alpha(2*m):alpha(2*m+1)-1)=-m;
            % It's all neatly ordered here, downgoing within every order
            %IMTAP(alpha(2*m):alpha(2*m+1)-1)=1:length(Vp);
        end
      % Duplicate for the positive order in case the region is axisymmetric
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
varns={G,V,EL,EM};
varargout=varns(1:nargout);

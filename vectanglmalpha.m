function [G,V]=vectanglmalpha(TH,L,srt)
% [G,V]=vectanglmalpha(TH,L,srt)
%
% Construction of tangent vector Slepian functions from the the Blm, and Clm
% vector spherical harmonics for spherical caps
%
% INPUT:
%
% TH    polar cap opening angle(s) ( no named regions yet)
% L     maximum spherical harmonic degree (no bandwidth at this point)
% srt   Should the Slepian functions be sorted? [default = 1 = yes]
%
% OUTPUT:
%
% G     Matrix containing the slepian coefficients for linear combinations
%       of the Elm and Flm. The first (Lmax+1)^2 coefficients are for the
%       Plm, the second (Lmax+1)^2-1 coefficients are for the Blm, and the
%       last (Lmax+1)^2-1 coefficients are for the Clm.
% V     concentration (eigen) values
%
% Last modified by plattner-at-alumni.ethz.ch, 05/05/2017

L=max(L);
defval('srt',1)

fname=fullfile(getenv('IFILES'),'VECTANGLMALPHA',...
		     sprintf('vectanglmalpha-%g-%g-%i-%i.mat',...
             max(TH),min(TH),max(L),min(L)));         
         
if exist(fname,'file')==2
    load(fname)
    disp(sprintf('Loading %s',fname))
else
    % For axisymmetric regions       
  
    mvec=0:L;
    sizesBC=max(L+1-max(mvec,1), zeros(size(mvec)) );
    siztot=2*sizesBC;
    
    deM=addmout(L);
    
    alpha=cumsum([1 siztot(1) gamini(siztot(2:end),2) ]);       
    
    % The redistribution is not even as complicated as it seems. The alpha
    % intervalls need to be exactly the size of the number of eigenvectors.
    % The right location with the ms is taken care of by EM==m 
    % (might need to make this faster?)

    % To simplify things and avoid mistakes: Treat the Blm and Clm 
    % coefficients the same way as the Elm coefficients and just remove the L=0 part if
    % necessary.
    
    % Initialize matrices  

    % Speeding up by preallocating the number of elements in the sparse
    % matrices
    % There are a bit less than sum_{m=0}^L (L+1-m)^2 elements.
    nelems=(L+1)^3 - L*(L+1)^2 +L*(L+1)*(2*L+1)/6;
    %GB=sparse((L+1)^2,2*(L+1)^2 - 2);
    GB=sparse([],[],[],(L+1)^2,2*(L+1)^2 - 2, nelems);
    % Treat Blm,Clm like Plm and remove the zero part later
    %GC=sparse((L+1)^2,2*(L+1)^2 - 2);
    GC=sparse([],[],[],(L+1)^2,2*(L+1)^2 - 2, nelems);

    V=zeros(1,2*(L+1)^2-2);
    disp('Calculating in parallel mode')
    
    % Now the same as in inoutgradvecglmalphaup: Calculate the individual
    % solutions for the m and put them back in the right place
    parfor mm=1:L+1         
            m=mm-1;      
            [~,~,~,Cp1,Vpp1]=capvectorslepian(L,TH,m,[],[],[],[],[],[],0);
            [~,~,~,Cp2,Vpp2]=capvectorslepian(L,TH,-m,[],[],[],[],[],[],0);
            %[Vpp,Cp]=inoutgradvecsdwcap(TH,Lin,Lout,m);
            Vp1{mm}=Vpp1;
            Vp2{mm}=Vpp2;
            sizBC=max(L+1-max(m,1),zeros(size(L+1-max(m,1))));         
            CB1{mm}=Cp1(1:sizBC,:);            
            CC1{mm}=Cp1(sizBC+1:end,:);
            CB2{mm}=Cp2(1:sizBC,:);            
            CC2{mm}=Cp2(sizBC+1:end,:);
    end 
    
    
    % To make distribution a bit simpler, add the L=0 row to CB and CC.
    % This is only necessary for m=0.
    CB1{1}=[zeros(1,size(CB1{1},2)); CB1{1}];
    CC1{1}=[zeros(1,size(CC1{1},2)); CC1{1}];
    CB2{1}=[zeros(1,size(CB2{1},2)); CB2{1}];
    CC2{1}=[zeros(1,size(CC2{1},2)); CC2{1}];
    
    % Distribute this at the right point in the huge matrix
    for m=0:max(L)   
        if m>0
            % Here you supply the negative orders
            GB(deM==-m,alpha(2*m):alpha(2*m+1)-1)=CB2{m+1};  
            GC(deM==-m,alpha(2*m):alpha(2*m+1)-1)=CC2{m+1};  
            V(alpha(2*m):alpha(2*m+1)-1)=Vp2{m+1};
        end
        % Duplicate for the positive order in case the region is axisymmetric
            GB(deM==m,alpha(2*m+1):alpha(2*m+2)-1)=CB1{m+1};     
            GC(deM==m,alpha(2*m+1):alpha(2*m+2)-1)=CC1{m+1};  
            V(alpha(2*m+1):alpha(2*m+2)-1)=Vp1{m+1};   
    end           

    GB=GB(2:end,:);
    GC=GC(2:end,:);
    G=[GB;GC];

    if srt
        [V,isrt]=sort(V,'descend');           
        G=G(:,isrt);
    end
    
    if exist('octave_config_info')
    	% Octave
    	save(fname,'G','V')
    else
    	% Matlab
    	save(fname,'G','V','-v7.3')
    end
     
    
end % end of calculation if not yet available
    
% Provide output
varns={G,V};
varargout=varns(1:nargout);       
    
    
    
    
    
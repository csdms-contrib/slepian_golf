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
% Last modified by plattner-at-alumni.ethz.ch, 05/09/2017

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

    V=nan(1,2*(L+1)^2-2);
    %disp('Calculating in parallel mode')
    
    % Now the same as in inoutgradvecglmalphaup: Calculate the individual
    % solutions for the m and put them back in the right place
    parfor mm=1:L+1  
    %for mm=1:L+1
            m=mm-1;      
            %[~,~,~,Cp1,Vpp1]=capvectorslepian(L,TH,m,[],[],[],[],[],[],0);
            %[~,~,~,Cp2,Vpp2]=capvectorslepian(L,TH,-m,[],[],[],[],[],[],0);            
            %[Vpp1,Cp1{mm},Vpp2,Cp2{mm}]=vectansdwcap(TH,L,m);
            [Vpppos,Cppos,Vppneg,Cpneg]=vectansdwcap(TH,L,m);
            %[Vpp,Cp]=inoutgradvecsdwcap(TH,Lin,Lout,m);
            Vppos{mm}=Vpppos;
            Vpneg{mm}=Vppneg;
            sizBC=max(L+1-max(m,1),zeros(size(L+1-max(m,1))));
            CBpos{mm}=Cppos(1:sizBC,:);            
            CCpos{mm}=Cppos(sizBC+1:end,:);
            CBneg{mm}=Cpneg(1:sizBC,:);            
            CCneg{mm}=Cpneg(sizBC+1:end,:);
    end    
    % To make distribution a bit simpler, add the L=0 row to CB and CC.
    % This is only necessary for m=0.
    CBpos{1}=[zeros(1,size(CBpos{1},2)); CBpos{1}];
    CCpos{1}=[zeros(1,size(CCpos{1},2)); CCpos{1}];
    CBneg{1}=[zeros(1,size(CBneg{1},2)); CBneg{1}];
    CCneg{1}=[zeros(1,size(CCneg{1},2)); CCneg{1}];    
    % Distribute this at the right point in the huge matrix   
    for m=0:L   
        % This is where the positive-m-blm and negative-m-clm need to
        % get switched. 
        % I think it is because in Plattner et al 2014 (ACHA), the
        % D_(lm,l'm') are only non-zero when m'=-m. at the same time,
        % B_lm,l'm=B_l-m,l'-m and C_lm,l'm=C_l-m,l'-m. So ultimately we are
        % solving for the -m-component of Blm when we are solving for the
        if m>0
            % Here you supply the negative orders

        % m-component of Clm and vice-versa.
            GB(deM==-m,alpha(2*m):alpha(2*m+1)-1)=CBneg{m+1};  
            GC(deM==m,alpha(2*m):alpha(2*m+1)-1)=CCneg{m+1}; 
            V(alpha(2*m):alpha(2*m+1)-1)=Vpneg{m+1};
        end
        % Duplicate for the positive order in case the region is axisymmetric  
            GB(deM==m,alpha(2*m+1):alpha(2*m+2)-1)=CBpos{m+1};    
            GC(deM==-m,alpha(2*m+1):alpha(2*m+2)-1)=CCpos{m+1};    
            V(alpha(2*m+1):alpha(2*m+2)-1)=Vppos{m+1};             
    end              
    GB=GB(2:end,:);
    GC=GC(2:end,:);
    
    G=[GB;GC];
    
    if exist('octave_config_info')
    	% Octave
    	save(fname,'G','V')
    else
    	% Matlab
    	save(fname,'G','V','-v7.3')
    end
     
    
end % end of calculation, if not yet available
    
if srt
    [V,isrt]=sort(V,'descend');           
    G=G(:,isrt);
end

% Provide output
varns={G,V};
varargout=varns(1:nargout);       
    
    
    
    
    
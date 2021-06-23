function [G,V]=glmalpharing(TH,L,srt)
% [G,V]=glmalpharing(TH,L,srt)
%
% Returns the scalar/radial Slepian coefficients for spherical caps or 
% rings.
%
% INPUT:
%
% TH    Semi-opening angle of the spherical cap or two semi-opening angles
%       for the ring between them
% L     Maximum spherical-harmonic degree
% srt   Would you like the output to be sorted based on eigenvalues
%       default=1=yes
%
% OUTPUT:
%
% G     Matrix, who's columns are the eigenvectors = spherical-harmonic
%       coefficients for the Slepian functions
% V     Vector of eigenvalues in the same order as the columns of G
%
% Last modified by plattner-at-alumni.ethz.ch, 5/11/2017

defval('srt',1)

fname=fullfile(getenv('IFILES'),'GLMALPHARING',...
      sprintf('glmalpharing-%f-%f-%i.mat',max(TH),min(TH),L));
      
      
if exist(fname,'file')==2
    load(fname)
    disp(sprintf('Loading %s',fname))      
else
    % We will need these to distribute the vectors calculated for each
    % m-block among the big matrix
    mvec=0:L;
    sizes=max(L+1-mvec, zeros(size(mvec)) );
    deM=addmout(L);    
    alpha=cumsum([1 sizes(1) gamini(sizes(2:end),2) ]);      
    
    % Initialize the sparse matrix
    nelems=(L+1)^3-L*(L+1)^2+L*(L+1)*(2*L+1)/6;
    G=sparse([],[],[],(L+1)^2,(L+1)^2,nelems);
    
    V=nan(1,(L+1)^2);
    % Precalculate all the m-blocks in a parallel for loop
    parfor mm=1:L+1
        m=mm-1;
        [Vm{mm},Cm{mm}]=sdwcapring(TH,L,m);
    end
        
    % Now assemble
    for m=0:L
        if m>0
            % Put the negative orders in first
            G(deM==-m,alpha(2*m):alpha(2*m+1)-1)=Cm{m+1};
            V(alpha(2*m):alpha(2*m+1)-1)=Vm{m+1};
        end
        G(deM==m,alpha(2*m+1):alpha(2*m+2)-1)=Cm{m+1};
        V(alpha(2*m+1):alpha(2*m+2)-1)=Vm{m+1};
    end
           
    try
    	% Matlab
    	save(fname,'G','V','-v7.3')
    catch
      % Octave
    	save(fname,'G','V')
    end
    
end


if srt
    [V,isrt]=sort(V,'descend');           
    G=G(:,isrt);
end

% Provide output
varns={G,V};
varargout=varns(1:nargout);       
    


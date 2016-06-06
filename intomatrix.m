function varargout=intomatrix(coeff,L,lzero)
% [coeffmatrix,coeff]=INTOMATRIX(coeff,L,lzero)
%
% Transforms a list of coefficients into the matrix m l, which can be 
% nicely plotted using IMAGEFNAN
% 
% INPUT: 
%
% coeff         The coefficients that should be transformed into the
%               matrix. Either as coefficient list or as [l m cos sin]
% L             The maximum degree of the coefficients
% lzero         does coeff contain l=0? (For example the radial case)
% 
% OUTPUT:
%
% coeffmatrix   The coefficients arranged in a matrix that is easy to plot
% coeff         The list of coefficients, in case the coefficients were
%               given as [l m cos sin] 
%
% See also VECTORSPECTRAL, IMAGEFNAN
%
% Last modified by plattner-at-alumni.ethz.ch, 02/27/2012


defval('lzero',1)

[demsz,delsz,mz,lmc,mzin,mzo]=addmon(L);

% Test if it is in the lmcosi format
if size(coeff,2)==4&&size(coeff,2)>1
    lmcosi=coeff;
    if ~lzero
        lmcosi=[nan(1,4);lmcosi];
    end
    % Transform into list    
    try
        lm=reshape(lmcosi(:,3:4),1,2*length(demsz));
        coeff=lm(mzo);   
    catch
        error(...
'Probably did not check the lzero option when applying to the tangential component')
    end
end

coeffmatrix=nan(L+1,2*(L+1)-1);
count=1;
for l=0:L
    % m=0
    coeffmatrix(l+1,L+1)=coeff(count);
    count=count+1;
    for m=1:l
        % The -m
        coeffmatrix(l+1,L+1-m)=coeff(count);
        count=count+1;
        % The +m
        coeffmatrix(l+1,L+1+m)=coeff(count);
        count=count+1;
    end
end

varns={coeffmatrix,coeff};
varargout=varns(1:nargout);

function [slepcoef,dataweights,dampweights]=itweighres(M,data,slepcoef,niter,lambda,slepcoef0)
% Iteratively reweighted residual solver based on 
% Farquharson and Oldenburg (1998)
%
% Use eq (23--26) in their paper. 
%
% 4/23/2018: Now with damping

defval('lambda',0)
defval('slepcoef0',zeros(size(slepcoef)))

defval('rthresh',(1e-4)*median(abs(data(:))))
defval('rthresh_s',(1e-4)*median(abs(slepcoef(:))))
% Here is the weighted residual iteration
weights=[];
weights_s=[];
for iter=1:niter       
    residuals=data-M'*slepcoef;      
    residuals(abs(residuals)<rthresh)=rthresh;

    % For the damping, we need different residuals
    residuals_s=slepcoef-slepcoef0;
    residuals_s(abs(residuals_s)<rthresh_s)=rthresh_s;

    % Put the inverse of the thresheld residuals as diagonal elements
    % in the weighting matrix    
    weights=1./abs(residuals);

    weights_s=1./abs(residuals_s);

    % Here I am reweighing the data points
    % I think I can't get around having a second matrix and a second
    % data vector, because I don't want to keep multiplying them
    % Here the matrix      
    Mp=M.*repmat(weights',size(M,1),1);   
    % Instead of creating the large matrix and then multiplying,
    % we will use a for loop. Matlab's JIT compilation will make it 
    % better and it will save a lot of memory
    %Mp=sparse(size(M,1),size(M,2));
    %%%Mp=NaN(size(M));
    %%%for i=1:size(M,1)        
    %%%    Mp(i,:)=M(i,:).*(weights');        
    %%%end
    % Here the right hand side
    datap=weights.*data;

    % In the iterative reweighted residual, the weights do not
    % accumulate. Maybe there is something like a cumulative iterative
    % reweighting, but that's not what we want to do here. 
    % % Redefine the new matrix and the new data from reweighting 
    % Wcum=W*Wcum; data=W*data; M=W*M;

    % And now solve the reweighted problem: 
    % Farquharson and Oldenburg (1998) eq. (23):
    % (G'*R*G)m=G'*R*d 
    % slepcoef = (M*Mp')\(M*datap);
    % With damping:
    % (lambda*Rs + G'*R*G)m = lambda*Rs*m0 + G'*R*d
    slepcoef = (lambda*diag(weights_s) + M*Mp')\(lambda*weights_s.*slepcoef0 + M*datap);
    fprintf('iteration step %d: weighted thresheld residual norm is %g\n',iter,norm(residuals)/sqrt(length(data)))
end
dataweights=weights;
dampweights=weights_s;

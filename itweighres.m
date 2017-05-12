function [slepcoef,dataweights]=itweighres(M,data,slepcoef,niter)
% Iteratively reweighted residual solver based on 
% Farquharson and Oldenburg (1998)

defval('rthresh',(1e-4)*median(abs(data(:))))
% Here is the weighted residual iteration
weights=[];
for iter=1:niter       
    residuals=data-M'*slepcoef;      
    residuals(abs(residuals)<rthresh)=rthresh;
    % Put the inverse of the thresheld residuals as diagonal elements
    % in the weighting matrix    
    weights=1./abs(residuals);

    % Here I am reweighing the data points
    % I think I can't get around having a second matrix and a second
    % data vector, because I don't want to keep multiplying them
    % Here the matrix      
    Mp=M.*repmat(weights',size(M,1),1);                
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
    slepcoef = (M*Mp')\(M*datap);
    fprintf('iteration step %d: weighted thresheld residual norm is %g\n',iter,norm(residuals)/sqrt(length(data)))
end
dataweights=weights;
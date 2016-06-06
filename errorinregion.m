function varargout=errorinregion(dat,dattrue,dom,c11cmn)
% [diffmat,relmat]=errorinregion(dat,dattrue,dom,c11cmn)
%
% Calculates the difference of two data sets within a specified region
%
% INPUT:
%
% dat       Calculated or reconstructed data
% dattrue   Known or true data set
% dom       Region inside which we want to know the error
% c11cmn    coordinates in which the data sets are given 
%           [default: 0 90 360 -90]
%
% OUTPUT:
%
% diffmat   Matrix with the same size as dat. Inside the region it contains
%           the difference between the data sets. Outside the region it is 
%           zero
% relmat    similar to diffmat but the difference divided by the truedat 
%           value at each point
%
% EXAMPLE:
% 
% errorinregion('demo1') plot the value 1 only inside Australia
%                                 and 0 everywhere else 
%
% errorinregion('demo2') plot the value 1 only inside Africa
%                                 and 0 everywhere else 
%
% See also VECTORSPECTRAL
%
% Last modified by plattner-at-alumni.ethz.ch, 02/28/2012

defval('c11cmn',[0 90 360 -90]);
defval('dattrue',zeros(size(dat)));

if ~isstr(dat)

XY=eval(sprintf('%s(10)',dom));
    
thNdom=90-max(XY(:,2)); 
thSdom=90-min(XY(:,2));

thNtot=90-c11cmn(4);
thStot=90-c11cmn(2);
thtot=linspace(thStot,thNtot,size(dat,1));

phWtot=c11cmn(1);
phEtot=c11cmn(3);
phtot=linspace(phWtot,phEtot,size(dat,2));

diffmat=zeros(length(thtot),length(phtot),size(dat,3));
relmat=zeros(length(thtot),length(phtot),size(dat,3));

% Now find the phi intervals for each theta in c11cmn
phint=dphregion(thtot,[],dom);

for i=1:length(thtot)
    % find, which data point belongs to the intervall over the continent
    for j=1:(size(phint,2)/2)        
        startpoint=max(find(phtot<=phint(i,2*j-1)))+1;
        endpoint=min(find(phtot>=phint(i,2*j)))-1;
        for k=startpoint:endpoint         
            for p=1:size(dat,3)
                diffmat(i,k,p)=dattrue(i,k,p)-dat(i,k,p);
                relmat(i,k,p)=(dattrue(i,k,p)-dat(i,k,p))/dattrue(i,k,p);
            end
        end     
    end
end

varns={diffmat,relmat};
varargout=varns(1:nargout);


elseif strcmp(dat,'demo1')    
    dom='australia'
    fthph1=ones(1810,3610);
    fthph2=zeros(1810,3610);
    thrange=[-90 90];
    phrange=[0 360];
    c11cmn=[phrange(1) thrange(2) phrange(2) thrange(1)];
    err=errorinregion(fthph1,fthph2,dom,c11cmn);
    imagefnan([phrange(1) thrange(2)],[phrange(2) thrange(1)],...
        err,'kelicol',[0 1],[],0.1,100); 
    hold on
    XY=eval(sprintf('%s(10)',dom));
    plot(XY(:,1),XY(:,2),'k')
    hold off
    
elseif strcmp(dat,'demo2')    
    dom='africa'    
    fthph1=ones(1810,3610);
    fthph2=zeros(1811,3610);
    thrange=[-90 90];
    phrange=[180 540];
    c11cmn=[phrange(1) thrange(2) phrange(2) thrange(1)];
    err=errorinregion(fthph1,fthph2,dom,c11cmn);   
    imagefnan([phrange(1) thrange(2)],[phrange(2) thrange(1)],...
        err,'kelicol',[0 1],[],0.1,100);  
    hold on
    XY=eval(sprintf('%s(10)',dom));
    plot(XY(:,1),XY(:,2),'k')
    hold off
    
end  
   
function varargout=kerneltancapm(TH,Lmax,m,method)
% [Mm,Mmm,Bm,Dm]=kerneltancapm(TH,Lmax,m,method);
% 
% Kernel for [Blm Cl-m]. Returns the matrix Mm=[Bm Dm;Dm Bm], where
% Bm=integral of Blm*Bl'm = integral of Cl-m*Cl'-m and
% Dm=integral of Blm*Cl'-m = integral of Cl'-m*Blm.
%
% Can also return the matrix for -m (Mmm)
%
% INPUT:
%
% TH        Opening angle of the polar cap or a range for "belt"-region, in
%           degrees
% Lmax      Bandwidth
% m         Specific order
% method    how to calculate the Bm. An interger ngl for Gauss-Legendre 
%           with ngl Gauss-Legendre points or 'paul' for Paul-Gaunt
%
% OUTPUT:
%
% Mm        Localization kernel of the tangential plane for spherical polar
%           cap for order m and all degrees max(1,m) to Lmax. Takes non-
%           orthogonality of Blm with Cl-m into account
% Mmm       Localization kernel for -m
% Bm        BlmBlm  = Cl-mCl-m component of the localization kernel
% Dm        BlmCl-m =-Bl-mClm component of the localization kernel (note
%           that Cl-mBlm = BlmCl-m)
%
%
% See also KERNELBM, CAPVECTORSLEPIAN, KERNELB, KERNELBP
%
% Last modified by plattner-at-alumni.ethz.ch, 05/17/2012

defval('TH',30)
defval('Lmax',18)
defval('m',0)
defval('method',max(2*Lmax+1,200))

if ~isstr(TH)    

% Also allow interval for TH
if length(TH)==1
    TH=[TH 0];
end
% Here we go with the formulation using x=cos(theta)to make it the same 
% as in kernelc. The north pole is th=0, the south pole is th=180
thS=TH(1);
thN=TH(2);%0;
x0=cos(thS*pi/180);
x1=cos(thN*pi/180);

% The l=0 component is purely radial, we hence start with l=1;
if Lmax==0
    Mm=[];
    Mmm=[];
else        
Lmin=max(1,abs(m));
len=Lmax-Lmin+1;
Dm=zeros(len);

% The Dm can be calculated analytically
% It is only nonzero if m~=0
if m~=0
for ind1=1:len
    for ind2=ind1:len                
        l1=Lmin+ind1-1;
        l2=Lmin+ind2-1;
        XlmTHl1=libbrecht(l1,[x0 x1],'sch',[],abs(m))*sqrt(2*l1+1)...
            /sqrt(2-(m==0));
        XlmTHl2=libbrecht(l2,[x0 x1],'sch',[],abs(m))*sqrt(2*l2+1)...
            /sqrt(2-(m==0));
        Dm(ind1,ind2)=-2*pi/sqrt(l1*(l1+1)*l2*(l2+1)).*...
            m*(XlmTHl1(1)*XlmTHl2(1) - XlmTHl1(2)*XlmTHl2(2));
    end
end
end

Dm=Dm+Dm'-diag(diag(Dm));
% To make this exactly equivalent to Tony's \ylm, i.e. undo what we
% did above here, taking the output of YLM and multiplying
Dm=Dm/4/pi;   

% Now get the Bm from kernelbm
Bm=kernelbm(TH,Lmax,abs(m),method);
Mm=[Bm Dm;Dm Bm];
Mmm=[Bm -Dm;-Dm Bm];

varns={Mm,Mmm,Bm,Dm};
varargout=varns(1:nargout);

end

elseif strcmp(TH,'demo1')
    % Calculate Dm using numerical integration and compare to the analytic
    % solution
    TH=30;
    Lmax=ceil(rand*30)
    m=ceil(rand*Lmax)*  (ceil(rand*2)-1.5)*2
    ngl=200;
    [~,~,~,Dm]=kerneltancapm(TH,Lmax,m);
    
    % Now the numerical integration
    Dmnum=zeros(size(Dm));
    % Preparing the Gauss-Legendre integration 
    thS=TH;
    thN=0;
    x0=cos(thS*pi/180);
    x1=cos(thN*pi/180);
    intv=[x0 x1];
    nGL=max(ngl,2*Lmax);
    [w,x,N]=gausslegendrecof(nGL,[],intv);
    disp(sprintf('%i Gauss-Legendre points and weights calculated',N))
    
    Lmin=max(1,abs(m));
    len=Lmax-Lmin+1;
    for ind1=1:len                
        l1=Lmin+ind1-1;
        [X dX]=libbrecht(l1,x(:)','sch',[],abs(m));  
         Xlm=( X*sqrt(2*l1+1)/sqrt(2-(m==0)))';
        dXlm=(dX*sqrt(2*l1+1)/sqrt(2-(m==0)))';
        
        for ind2=1:len                            
            l2=Lmin+ind2-1;   
            [Xp dXp]=libbrecht(l2,x(:)','sch',[],abs(m));  
             Xlpm=( Xp*sqrt(2*l2+1)/sqrt(2-(m==0)))';
            dXlpm=(dXp*sqrt(2*l2+1)/sqrt(2-(m==0)))';
            
            % The sin(acos(x)) transformation is necessary, because the
            % derivative is with respect to theta=cos(x).
            if(m>0)           
                Dmnum(ind1,ind2)=2*(-m)/sqrt(l1*(l1+1)*l2*(l2+1))*(...
                    w(:)'*(dXlm(:).* Xlpm(:)./sin(acos(x)))*pi+... 
                    w(:)'*( Xlm(:).*dXlpm(:)./sin(acos(x)))*pi ...
                );
            else
                Dmnum(ind1,ind2)=2*(-m)/sqrt(l1*(l1+1)*l2*(l2+1))*(...
                    w(:)'*(dXlm(:).* Xlpm(:)./sin(acos(x)))*pi+...
                    w(:)'*( Xlm(:).*dXlpm(:)./sin(acos(x)))*pi...
                );
            end
        end
    end
    
    % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
    % did above here, taking the output of YLM and multiplying
    Dmnum=Dmnum/4/pi;   

    disp(sprintf(...
'Average relative difference between analytical and numerical D kernel is %g'...
,norm(Dmnum-Dm)./norm(Dm)));
   
end

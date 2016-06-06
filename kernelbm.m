function varargout=kernelbm(TH,L,m,method)
% Bm=kernelbm(TH,Lmax,m,method)
%  
% integral Blm*Bl'm = integral Clm*Cl'm for a SINGLE order m and all 
% degrees l,l' between m and the maximum degree Lmax.
% This is NOT the kernel for a single order. Each order m interacts with
% -m. Therefore, use KERNELTANCAPM to set up the kernel for a single order.
%
% INPUT:
% 
% TH        Angular extent of the spherical cap or a range for "belt"-
%           region, in degrees
% L         Bandwidth
% m         Angular order
% method    an interger ngl for Gauss-Legendre with ngl Gauss-Legendre 
%           points or 'paul' for Paul-Gaunt
%
% OUTPUT:
%
% Bm         Localization kernel for order m
%
% EXAMPLE:
%
% kernelbm('demo1') Compares Bm for different numbers of Gauss points to
%                   the assembly using the method of Paul
%
% See also KERNELTANCAPM, CAPVECTORSLEPIAN, KERNELB, KERNELBP
%
% Last modified by plattner-at-alumni.ethz.ch, 04/22/2015

defval('TH',30)
defval('Lmax',18)
defval('m',0)
defval('method',max(200,2*max(L)+1))

bp=length(L)==2;% new, amp 4/22/2015
lmin=max(m,bp*min(L)); % new, amp 4/22/2015
Lmax=max(L); % new, amp 4/22/2015

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

% B_{l1,m,l2,m} = B_{l1,-m,l2,-m}, therefore working with the absolute 
% value of m 
m=abs(m);
% The l=0 component is purely radial, we hence start with l=1;
if Lmax==0
    Bm=[];
else    



Lmin=max(1,m);
Lmin=max(lmin,Lmin); % new, amp 4/22/2015
len=Lmax-Lmin+1;
Bm=zeros(len); 

if ~isstr(method)
    ngl=method;    
   
    % Approximate Nyquist degree is $\pi/d\theta$
    if ngl < (Lmax+1) & ngl~=0
        error('Sample finer to avoid aliasing')
    end

    tic
    [w,x,N]=gausslegendrecof(ngl,[],[x0 x1]);
    
    % Calculate the Xlm and the dXlm for L=m to Lmax and for all 
    % Gauss-Legendre points    
     Xlm=repmat(NaN,length(x),len);
    dXlm=repmat(NaN,length(x),len);
    ind=0;       
    for ind=1:len 
        l=Lmin+ind-1;
        [X dX]=libbrecht(l,x(:)','sch',[],m);   
        Xlm(:,ind)=(X*sqrt(2*l+1)/sqrt(2-(m==0)))';
        dXlm(:,ind)=(dX*sqrt(2*l+1)/sqrt(2-(m==0)))';         
    end
           
    % Calculate the Legendre products for all combinations of l1 and l2
    BGL=repmat(NaN,length(x),((len^2)+len)/2);
    index=0;
    
    %h=waitbar(0,'KERNELBM: Calculating all Legendre products');
    for ind1=1:len
        for ind2=ind1:len
            index=index+1;            
            l1=Lmin+ind1-1;
            l2=Lmin+ind2-1;
            BGL(:,index)=2*pi/sqrt(l1*(l1+1)*l2*(l2+1))*(...
                dXlm(:,ind1).*dXlm(:,ind2) + ...
                (m*Xlm(:,ind1)./sin(acos(x))).*...
                (m*Xlm(:,ind2)./sin(acos(x))) );                    
            %waitbar(index/(len*(len+1)/2),h)
        end
    end
    %delete(h)
    
    Bvec=BGL'*w;
    
    % Transform the vector into an upper triangular matrix
    index=0;
    for l1=1:len
        for l2=l1:len
        index=index+1;
        Bm(l1,l2)=Bvec(index);
        end
    end
    Bm=Bm+Bm'-diag(diag(Bm));
    
    toc
    
    % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
    % did above here, taking the output of YLM and multiplying
    Bm=Bm/4/pi;                         

elseif strcmp(method,'paul')  
    % Using the method of Paul-Gaunt.           
    tic                 
    % Now assembling all Integrals of a single associated Legendre
    % function at the theta intervals. 2*Lmax is needed because of the
    % formula of Gaunt (need the correctly normalized ones). 
    % Allow ranges for TH and ultimately for x
    [~,Itab]=paul(2*Lmax,[x0 x1]);
    Itab=Itab(:,1)-Itab(:,2);
    
    % Next, load the wigner symbols up to level 2*Lmax
    try
     [~,C0,S0,LLoad0]=zeroj(0,0,2*Lmax);
    catch
     wignercycle(2*Lmax,0,0);
     [~,C0,S0,LLoad0]=zeroj(0,0,2*Lmax);
    end

    % The part here is repeated mutatis mutandis from threej. 
    % It is done to load the smallest possible database.
    % First check which databases are available:
    try
       Els=ls2cell(fullfile(getenv('IFILES'),...
                       'WIGNER','WIGNER3J-*-sym.mat'));
    catch
       Els=[];
       EL=[];
    end
    for index=1:length(Els)
       EL(index)=str2num(rindeks(parse(Els{index},'-'),2));
    end
    EL=sort(EL);
    % Bandwidth of the database; keep this as low as possible
    % Need to provide for empties if not found
    fmax=find(2*Lmax<=EL);
    if ~isempty(fmax)
        LLoad3=EL(indeks(fmax,1));
        % Check, identify and load database
        if 2*Lmax>LLoad3
            disp('Generating wigner3j database ...')
            wignercycle(2*Lmax,1,1);
            disp('... done.')                  
        end
        wign=sprintf('%s/WIGNER3J-%i-sym',...
                fullfile(getenv('IFILES'),'WIGNER'),LLoad3);
        disp(sprintf('Loading %s',wign))
        load(wign)
    else
        disp('Generating wigner3j database ...')
        wignercycle(2*Lmax,1,1);
        disp('... done.')
        wign=sprintf('%s/WIGNER3J-%i-sym',...
                    fullfile(getenv('IFILES'),'WIGNER'),2*Lmax);
        load(wign)
    end
    
    % Loading complete, now assemble the matrix
    index=0;    
    h=waitbar(0,'KERNELBM: Calculating all Legendre products');
    for ind1=1:len
    for ind2=ind1:len
        index=index+1;            
        l1=Lmin+ind1-1;
        l2=Lmin+ind2-1;

        % The Ilk coefficients for dX
        a1_l1m1=-sqrt((l1+m).*(l1-m+1))/2;
        a2_l1m1= sqrt((l1-m).*(l1+m+1))/2;
        a1_l2m2=-sqrt((l2+m).*(l2-m+1))/2;
        a2_l2m2= sqrt((l2-m).*(l2+m+1))/2;

        % The Ilk coefficients for mX
        b1_l1m1=-sqrt((2*l1+1)/(2*l1-1))*sqrt((l1+m).*(l1+m-1))/2;
        b2_l1m1=-sqrt((2*l1+1)/(2*l1-1))*sqrt((l1-m).*(l1-m-1))/2;
        b1_l2m2=-sqrt((2*l2+1)/(2*l2-1))*sqrt((l2+m).*(l2+m-1))/2;
        b2_l2m2=-sqrt((2*l2+1)/(2*l2-1))*sqrt((l2-m).*(l2-m-1))/2;

        IdXdX=... 
         (a1_l1m1*a1_l2m2*...
            calcGaunt(l1,m-1,l2,m-1,LLoad0,C0,S0,w3js,Itab)+...                                
          a1_l1m1*a2_l2m2*...
            calcGaunt(l1,m-1,l2,m+1,LLoad0,C0,S0,w3js,Itab)+...
          a2_l1m1*a1_l2m2*...
              calcGaunt(l1,m+1,l2,m-1,LLoad0,C0,S0,w3js,Itab)+...
          a2_l1m1*a2_l2m2*...
              calcGaunt(l1,m+1,l2,m+1,LLoad0,C0,S0,w3js,Itab) ...
         )*sqrt(2*l1+1)*sqrt(2*l2+1)*sqrt(2-(m==0))*sqrt(2-(m==0));


        ImXmX=...
         (b1_l1m1*b1_l2m2*...
              calcGaunt(l1-1,m-1,l2-1,m-1,LLoad0,C0,S0,w3js,Itab)+...                    
          b1_l1m1*b2_l2m2*...
              calcGaunt(l1-1,m-1,l2-1,m+1,LLoad0,C0,S0,w3js,Itab)+...                    
          b2_l1m1*b1_l2m2*...
              calcGaunt(l1-1,m+1,l2-1,m-1,LLoad0,C0,S0,w3js,Itab)+...                     
          b2_l1m1*b2_l2m2*...
              calcGaunt(l1-1,m+1,l2-1,m+1,LLoad0,C0,S0,w3js,Itab) ...                     
         )*sqrt(2*(l1-1)+1)*sqrt(2*(l2-1)+1)*sqrt(2-(m==0))*sqrt(2-(m==0));                   

        waitbar(2*index/((len^2)+len),h);
        Bm(ind1,ind2)=2*pi*(IdXdX+ImXmX)/sqrt(l1*(l1+1)*l2*(l2+1))...
            /(2-(m==0));
        Bm(ind2,ind1)=Bm(ind1,ind2);
    end
    end
    delete(h)  
    toc    
    % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
    % did above here, taking the output of YLM and multiplying
    Bm=Bm/4/pi;                  
else
    error('Select valid method');    
end
end

varns={Bm};
varargout=varns(1:nargout);

% The demos
elseif strcmp(TH,'demo1')  
    % Calculate Bm with different methods and compare them
    TH=ceil(rand*90);
    Lmax=ceil(rand*(20));
    m=floor(rand*(Lmax+1));
    ngl=50:50:500;     
    Bmpaul=kernelbm(TH,Lmax,m,'paul');
    for i=1:length(ngl)
        Bmgl=kernelbm(TH,Lmax,m,ngl(i));        
        diff(i)=sum(sum(abs(Bmgl-Bmpaul)));
    end  
    plot(ngl,diff)
    xlabel('Number of Gauss points')
    ylabel('Error compared to solution using paul')
    
end

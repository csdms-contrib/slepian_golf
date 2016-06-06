function varargout=calcilk(L,x,pm,znorm,CP)
% [mXlm,dXlm,mXlmcompare,dXlmcompare]=calcilk(L,x,pm,znorm,CP)
%
% Uses Ilk Lemma 1 and 2 to calculate m*X_{lm}/sin from the X_{l-1}. And
% dX_{lm} from the X_{l}. On demand also calculates m*X/sin and dX the
% classical way to compare.
%
% INPUT:
%
% L             Maximum angular degree [default=10]
% x             the cos(theta) values to evaluate [default -0.99 to 1]
% pm            use (-1)^m factor in the Xlm? [default: 0]
% znorm         use sqrt(2) renormalization for m=0? [default 0]
% CP            precalculated  Xlm values, or linear functionals of the Xlm
%               in order to apply ilk directly to the functional values of 
%               the Xlm (for example integrals of the Xlm) 
%               [default: empty (i.e. don't use)]
%
% OUTPUT:
%
% mXlm          the m*X_{lm}/sin for the chosen degree L
% dXlm          the dX_{lm} for the chosen degree L
% mXlmcompare   the m*X_{lm}/sin calculated the classical way (is only 
%               calculated on demand)
% dXlmcompare   the dX_{lm} calculated the classical way (is only 
%               calculated on demand)
%
% EXAMPLE:
%
% calcilk('demo1') Compare dX and mX/sin calculated with Ilk and 
%                           in the classical way
%
% calcilk('demo2') Compare the integrals of dX and mX/sin
%                           calculated with Paul/Ilk and the classical way
%                           using Gauss-Legendre
%
% See also XYZ2BLMCLM
%
% Last modified by plattner-at-alumni.ethz.ch, 02/28/2012


defval('L',10)
defval('x',-0.99:0.01:1);
defval('pm',0)
defval('znorm',0)
defval('CP',[])
x=x(:)';
N=length(x);

if ~ischar(L)

mXlm=NaN(length(x),addmup(L));
mXlm(:,1)=zeros(length(x),1);
dXlm=NaN(length(x),addmup(L));
dXlm(:,1)=zeros(length(x),1);
mXlmcompare=[];
dXlmcompare=[];


if(nargout>2)
    mXlmcompare=NaN(length(x),addmup(L));
    mXlmcompare(:,1)=zeros(length(x),1);
    dXlmcompare=NaN(length(x),addmup(L));
    dXlmcompare(:,1)=zeros(length(x),1);
end

in1=1;
in2=3;
oldin1=0;
oldin2=1;
for l=1:L
    if isempty(CP)   
        Pm1=legendre(l-1,x,'sch')*sqrt(2*(l-1)+1);
        P=legendre(l,x,'sch')*sqrt(2*l+1);
    else
        Pm1=CP(oldin1+1:oldin2,:);
        P=CP(in1+1:in2,:);               
    end   
    % Normalizations
    % Pm1:
    nrm=ones(size(Pm1));
    nrm(1,:)=sqrt(2);
    Pm1=Pm1.*nrm;
    Pm1=[Pm1;zeros(2,size(Pm1,2))];
    Pm1=[(-1)*Pm1(2,:);Pm1];    
    Xm1=Pm1.*repmat((-1).^(-1:size(Pm1,1)-2),size(Pm1,2),1)';
    % P:
    nrm=ones(size(P));
    nrm(1,:)=sqrt(2);
    P=P.*nrm;
    P=[P;zeros(1,size(P,2))];
    P=[(-1)*P(2,:);P];    
    X=P.*repmat((-1).^(-1:size(P,1)-2),size(P,2),1)';    
    m=(0:l)';
    
    % Ilk Lemma 1
    a1=-repmat(sqrt((l+m).*(l-m+1))/2,1,N);
    a2= repmat(sqrt((l-m).*(l+m+1))/2,1,N);  
    dXlm(:,in1+1:in2)=(a1.*X(1:end-2,:)+a2.*X(3:end,:))';
    
    % Ilk Lemma 2
    b1=-repmat(sqrt((2*l+1)/(2*l-1))*sqrt((l+m).*(l+m-1))/2,1,N);
    b2=-repmat(sqrt((2*l+1)/(2*l-1))*sqrt((l-m).*(l-m-1))/2,1,N);  
    mXlm(:,in1+1:in2)=(b1.*Xm1(1:end-2,:)+b2.*Xm1(3:end,:))';
    
    if(nargout>2&&isempty(CP))
        Pcomp=libbrecht(l,x,'sch');
        Xcomp=Pcomp.*repmat((-1).^(0:size(Pcomp,1)-1),size(Pcomp,2),1)'...
            *sqrt(2*l+1);
        div_sinx=repmat(1./sin(acos(x)),length(m),1);
        mXlmcompare(:,in1+1:in2)=(repmat(m,1,size(Pcomp,2)).*...
            Xcomp(1:end,:).*div_sinx)';
        [~,dPcomp]=libbrecht(l,x,'sch');
        % For the dXlm comparison, the normalization of the m=0 case does
        % matter because we are not multiplying with m as in the Xlm case
        nrm=ones(size(dPcomp));
        nrm(1,:)=sqrt(2);
        dPcomp=dPcomp.*nrm;
        dXlmcompare(:,in1+1:in2)=(dPcomp').*...
           repmat((-1).^(0:size(dPcomp,1)-1),size(dPcomp,2),1)*sqrt(2*l+1);
    end
   
    if pm==0
        mXlm(:,in1+1:in2)=mXlm(:,in1+1:in2).*...
            repmat((-1).^(0:(in2-in1-1)),N,1);
        dXlm(:,in1+1:in2)=dXlm(:,in1+1:in2).*...
            repmat((-1).^(0:(in2-in1-1)),N,1);
        if (nargout>2&&isempty(CP))
            mXlmcompare(:,in1+1:in2)=mXlmcompare(:,in1+1:in2).*...
                repmat((-1).^(0:(in2-in1-1)),N,1);
            dXlmcompare(:,in1+1:in2)=dXlmcompare(:,in1+1:in2).*...
                repmat((-1).^(0:(in2-in1-1)),N,1);
        end
    end
    
    if znorm==0
        rnrm=ones(size(mXlm(:,in1+1:in2)));
        rnrm(:,1)=1/sqrt(2);
        mXlm(:,in1+1:in2)=mXlm(:,in1+1:in2).*rnrm;
        dXlm(:,in1+1:in2)=dXlm(:,in1+1:in2).*rnrm;
        if (nargout>2&&isempty(CP))
            mXlmcompare(:,in1+1:in2)=mXlmcompare(:,in1+1:in2).*rnrm;
            dXlmcompare(:,in1+1:in2)=dXlmcompare(:,in1+1:in2).*rnrm;
        end
        
    end
    
    oldin1=in1;
    oldin2=in2;    
    in1=in2;
    in2=in1+l+2;  
    
end

varns={mXlm,dXlm,mXlmcompare,dXlmcompare};
varargout=varns(1:nargout);

elseif strcmp(L,'demo1')
    L=round(rand*50);
    x=-0.99:0.01:1;
    [mXlm,dXlm,mXlmcompare,dXlmcompare]=calcilk(L,x);
    
    m=floor(rand*(L+1));
    
    subplot(2,2,1)
    plot(x,mXlm(:,m),x,mXlmcompare(:,m),'--r')
    legend('Ilk','Classic')
    title('mXlm')
    subplot(2,2,2)
    plot(x,mXlm(:,m)-mXlmcompare(:,m))
    title('Difference mXlm')
    
    subplot(2,2,3)
    plot(x,dXlm(:,m),x,dXlmcompare(:,m),'--r')
    legend('Ilk','Classic')
    title('dXlm')
    subplot(2,2,4)
    plot(x,dXlm(:,m)-dXlmcompare(:,m))
    title('Difference dXlm')
    
    
elseif strcmp(L,'demo2')
    L=floor(rand*31);
    x0=-0.99:0.1:1;
    m=floor(rand*(L+1));
    Itab=paul(L,x0);
    % Now calculate int_x0^1 mXlm(x)/sin(x) \,dx and int_x0^1 dXlm(x) \,dx
    % in two different ways: with ilk and with numerical integration
    % With ilk (all integrals for all L and m)
    [mXilk,dXilk]=calcilk(L,x0,[],[],Itab);
    % With numerical integration
 
    for i=1:length(x0)
        [wgl,xgl]=gausslegendrecof(100,[],[x0(i) 1]);
        [X,dX]=libbrecht(L,xgl,'sch');
        % The sqrt(2*L+1)/sqrt(2*(L-1)+1) factor is necessary because the
        % ilk calculation is written for the Xlm and here we have the Plm.
        % In dX this does not matter because the ilk linear combination
        % only uses Xlms of the same level (they are all multiplied with
        % the same 1/(2l+1). But in the mX case we take linear combinations
        % of one lower degree l-1.
        mXgl(i)=wgl(:)'*(m*X(m+1,:)'./sin(acos(xgl(:))))*sqrt(2*L+1)...
            /sqrt(2*(L-1)+1);
        dXgl(i)=wgl(:)'*(dX(m+1,:)');
    end
    
    entry=L*(L+1)/2+m+1;
    
    subplot(2,2,1)
    plot(x0,mXilk(:,entry),x0,mXgl,'--r')
    legend('Ilk+Paul','GL')
    title('mX')
    subplot(2,2,2)
    plot(x0,mXilk(:,entry)-mXgl')
    title('Difference mX')
    
    subplot(2,2,3)
    plot(x0,dXilk(:,entry),x0,dXgl,'--r')
    legend('Ilk+Paul','GL')
    title('dXlm')
    subplot(2,2,4)
    plot(x0,dXilk(:,entry)-dXgl')
    title('Difference dX')

    
end
    
    
    
    


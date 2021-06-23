function varargout=kernelb(Lmax,dom,pars,method,rotb,nstripes)
% [K,B,D,XY]=kernelb(Lmax,dom,pars,method,rotb,nstripes)
%  
% Vector spherical harmonic localization kernel for the tangential space
% for a all degrees l between 0 and Lmax and for all orders -m<=l<=m.
% Returns the full tangential Kernel and the B = C and D component.
%
% INPUT:
%
% Lmax      Bandwidth
% dom       Either the region 'england', 'eurasia',  'namerica' [default], 
%           'australia', 'greenland', 'africa', 'samerica', 'amazon', 
%           'gpsnamerica', 'antarctica', 'alloceans', 'orinoco'
%           OR: 'patch' spherical patch with specs in 'pars'
%           OR: 'sqpatch' square patch with [thN thS phW phE] in 'pars'
%           OR: [lon lat] an ordered list defining a closed curve in
%           degrees
% pars      [th0,ph0,thR] for 'patch'
%                th0  Colatitude of the cap center, in radians
%                ph0  Longitude of the cap center, in radians
%                thR  Radius of the cap, in radians
%           OR: [thN thS phW phE] for 'sqpatch'
%           OR: N  splining smoothness for geographical regions 
%               [default: 10]
%           OR: the string with the name you want the result saved as
% method    an interger for Gauss-Legendre with the indicated number of
%           Gauss-Legendre points or 'paul' for Paul-Gaunt
% rotb      0 That's it, you got it [default: 0]
%           1 For, e.g., 'antarctica', if you were given rotated 
%           coordinates to make the integration procedure work, this option
%           makes sure that the kernel matrix reflects this. If not, you 
%           have to apply counterrotation after diagonalizing in 
%           LOCALIZATION.
% nstripes  Number of equally wide longitudinal stripes to calculate the
%           integral over the continent in the 'paul' case
%
% OUTPUT:
%
% K         Localization kernel for the tangential plane
%           indexed as: degree  [1  1  1  2  2  2  2  2]
%                       order   [0 -1  1  0 -1  1 -2  2]
% B         The b*b=c*c part of the localization kernel 
% D         The b*c=-c*b part of the localization kernel
% XY        The outlines of the region into which you are localizing
%
% EXAMPLE:
%
% kernelb('demo1') Plot the eigenvalues for Africa
% 
% kernelb('demo2') Compare the results calculated with Gauss-Legendre to 
%                  the Paul results  
%
% See also VECTORSLEPIAN, BLMCLM2XYZ, KERNELBP, KERNELTANCAPM
%
% On 06/29/2017, plattner-at-alumni.ethz.ch made
% Antartica rotate back by default
%  
% Last modified by plattner-at-alumni.ethz.ch, 07/14/2017

defval('Lmax',18)
defval('dom','namerica')
defval('method',200)
defval('rotb',0)
defval('nstripes',200)

if strcmp(dom,'antarctica')
  rotb = 1;
end

if ~isstr(Lmax)
% Generic path name that I like
filoc=fullfile(getenv('IFILES'),'KERNELB');    

if ~(isstr(dom) || length(dom)>2 )
    error('Use CAPVECTORSLEPIAN to calculate Slepians for a polar cap')
else

    ngl=method;
    
    switch dom
        % If the domain is a square patch
    case 'sqpatch'
        fnpl=sprintf('%s/%s-%i-%i-%i-%i-%i-tang.mat',filoc,dom,Lmax,...
                round(pars(1)*180/pi),round(pars(2)*180/pi),...
                round(pars(3)*180/pi),round(pars(4)*180/pi));
        % If the domain is a spherical patch
    case 'patch'
        fnpl=sprintf('%s/%s-%i-%i-%i-%i-tang.mat',filoc,dom,Lmax,...
               round(pars(1)*180/pi),round(pars(2)*180/pi),...
               round(pars(3)*180/pi));
        % If the domain is a named region or a closed contour
    otherwise
        if ischar(dom)
	  h=dom;
	else 
	  try
	    h=hash(dom,'sha1');
    catch
      h=builtin('hash','sha1',dom);
	  end
	end
        fnpl=sprintf('%s/WREG-%s-%i-tang.mat',filoc,h,Lmax);
        % For some of the special regions it makes sense to distinguish
        % If it gets rotb=1 here, it doesn't in LOCALIZATION
        if strcmp(dom,'antarctica') && rotb==1 
            fnpl=sprintf('%s/WREG-%s-%i-%i-tang.mat',filoc,dom,Lmax,rotb);
        end
    end
    
    
    if exist(fnpl,'file')==2 && (~isstr(ngl) || strcmp(ngl,'paul'))
        load(fnpl)
        disp(sprintf('%s loaded by KERNELB',fnpl))
    else       
    % Calculate it    
    if strcmp(dom,'patch')
        defval('pars',[90 75 30]*pi/180);
        % For future reference 
        th0=pars(1); ph0=pars(2); thR=pars(3);
        if th0==0
        disp('Not meant for polar caps, running anyhow for comparison')
        %error('Not for polar caps! Use GRUNBAUM or SDWCAP instead')
        % BUT IN COMPARING, NOTE THAT THE SIGN MAY BE OFF
        end
        if thR>th0
        %disp('Really, should be putting in the GRUNBAUM call here')
        error('Not for near-polar caps! Use polar cap option, then rotate')
        end
        % Convert all angles to degrees for CAPLOC only
        [lon,lat]=caploc([ph0 pi/2-th0]/pi*180,thR/pi*180,100,1);
        % Northern and Southern points, in radians
        thN=(th0-thR);
        thS=(th0+thR);
        XY=[lon lat];
    elseif strcmp(dom,'sqpatch')
        defval('pars',[30 90 10 90]*pi/180);
        thN=pars(1); thS=pars(2); phW=pars(3); phE=pars(4);
        XY=[phW pi/2-thN ; phW pi/2-thS ; phE pi/2-thS ; ...
        phE pi/2-thN ; phW pi/2-thN]*180/pi;
    else
        if isstr(dom)
            % If it's a named geographical region or a coordinate boundary
            defval('pars',10);
            % Run the named function to return the coordinates
            if strcmp(dom,'antarctica') && rotb==1
                % Return the rotation parameters also, to undo later
                [XY,lonc,latc]=eval(sprintf('%s(%i)',dom,pars));
            else
                % Don't, the result will be the kernel for the rotated dom                
                XY=eval(sprintf('%s(%i)',dom,pars));
            end
	else
	  XY=dom;
        end
        thN=90-max(XY(:,2)); thN=thN*pi/180;
        thS=90-min(XY(:,2)); thS=thS*pi/180;
    end

    [dems,dels,mz,lmc,mzi,mzo,bigm,bigl]=addmon(Lmax);
    dimK=(Lmax+1)^2; lenm=length(dems);
    B=repmat(NaN,dimK,dimK);
    D=repmat(NaN,dimK,dimK);

    % Here comes the decision if it is GL or paul
    if isstr(ngl) && strcmp(ngl,'paul')
        % Using the method of Paul-Gaunt.           

        % choose an interval for the theta values (evaluation of the 
        % integrals over sinsin etc.). Here: equidistant           
        thvals=linspace(cos(thS),cos(thN),nstripes+1);                   

        % Now assembling all Integrals of a single associated Legendre
        % function at the theta intervals. 2*Lmax is needed because of the
        % formula of Gaunt (need the correctly normalized ones).
        [~,Itab]=paul(2*Lmax,thvals);   

        % Next, load the wigner symbols up to level 2*Lmax
        try
         [~,C0,S0,LLoad0]=zeroj(0,0,2*Lmax);
        catch
         wignercycle(2*Lmax,0,0);
         [~,C0,S0,LLoad0]=zeroj(0,0,2*Lmax);
        end

        % The part here is repeated mutatis mutandis from threej. It is 
        % done to load the smallest possible database.
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
            %defval('L',EL(indeks(fmax,1)))
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
       
        IdXdX=repmat(NaN,nstripes+1,((lenm^2)+lenm)/2);
        ImXmX=repmat(NaN,nstripes+1,((lenm^2)+lenm)/2);
        IdXmX=repmat(NaN,nstripes+1,((lenm^2)+lenm)/2);
        ImXdX=repmat(NaN,nstripes+1,((lenm^2)+lenm)/2);
        dubs=repmat(2,lenm,1); dubs(mz)=1; comdi=[];
        index=0;
        h=waitbar(0,...
            'KERNELB: Calculating product integrals using Paul-Gaunt');
        for lm1dex=1:lenm
           l1=dels(lm1dex);
           m1=dems(lm1dex);
           pos1=1/2*(l1)^2+1/2*l1+m1+1;
           comdi=[comdi ; dubs(lm1dex:lenm)];
           % Note that the last index will be ((lenm^2)+lenm)/2
           for lm2dex=lm1dex:lenm
            l2=dels(lm2dex);
            m2=dems(lm2dex);
            index=index+1;
            pos2=1/2*(l2)^2+1/2*l2+m2+1;            

            % The Ilk coefficients for dX
            a1_l1m1=-sqrt((l1+m1).*(l1-m1+1))/2;
            a2_l1m1= sqrt((l1-m1).*(l1+m1+1))/2;
            a1_l2m2=-sqrt((l2+m2).*(l2-m2+1))/2;
            a2_l2m2= sqrt((l2-m2).*(l2+m2+1))/2;

            % The Ilk coefficients for mX
            b1_l1m1=-sqrt((2*l1+1)/(2*l1-1))*sqrt((l1+m1).*(l1+m1-1))/2;
            b2_l1m1=-sqrt((2*l1+1)/(2*l1-1))*sqrt((l1-m1).*(l1-m1-1))/2;
            b1_l2m2=-sqrt((2*l2+1)/(2*l2-1))*sqrt((l2+m2).*(l2+m2-1))/2;
            b2_l2m2=-sqrt((2*l2+1)/(2*l2-1))*sqrt((l2-m2).*(l2-m2-1))/2;
                                   
            IdXdX(:,index)=... 
                (a1_l1m1*a1_l2m2*...
                    calcGaunt(l1,m1-1,l2,m2-1,LLoad0,C0,S0,w3js,Itab)+...                                
                 a1_l1m1*a2_l2m2*...
                    calcGaunt(l1,m1-1,l2,m2+1,LLoad0,C0,S0,w3js,Itab)+...
                 a2_l1m1*a1_l2m2*...
                    calcGaunt(l1,m1+1,l2,m2-1,LLoad0,C0,S0,w3js,Itab)+...
                 a2_l1m1*a2_l2m2*...
                    calcGaunt(l1,m1+1,l2,m2+1,LLoad0,C0,S0,w3js,Itab) ...
                )*sqrt(2*l1+1)*sqrt(2*l2+1)...
                    *sqrt(2-(m1==0))*sqrt(2-(m2==0));
             
            ImXmX(:,index)=...
                (b1_l1m1*b1_l2m2*...
                    calcGaunt(l1-1,m1-1,l2-1,m2-1,LLoad0,C0,S0,w3js,Itab)+...                    
                 b1_l1m1*b2_l2m2*...
                    calcGaunt(l1-1,m1-1,l2-1,m2+1,LLoad0,C0,S0,w3js,Itab)+...                    
                 b2_l1m1*b1_l2m2*...
                    calcGaunt(l1-1,m1+1,l2-1,m2-1,LLoad0,C0,S0,w3js,Itab)+...                     
                 b2_l1m1*b2_l2m2*...
                    calcGaunt(l1-1,m1+1,l2-1,m2+1,LLoad0,C0,S0,w3js,Itab) ...                     
                )*sqrt(2*(l1-1)+1)*sqrt(2*(l2-1)+1)...
                    *sqrt(2-(m1==0))*sqrt(2-(m2==0));                   
                              
            IdXmX(:,index)=...
                (a1_l1m1*b1_l2m2*...
                    calcGaunt(l1,m1-1,l2-1,m2-1,LLoad0,C0,S0,w3js,Itab)+...                     
                 a1_l1m1*b2_l2m2*...
                    calcGaunt(l1,m1-1,l2-1,m2+1,LLoad0,C0,S0,w3js,Itab)+...
                 a2_l1m1*b1_l2m2*...
                    calcGaunt(l1,m1+1,l2-1,m2-1,LLoad0,C0,S0,w3js,Itab)+...
                 a2_l1m1*b2_l2m2*...
                    calcGaunt(l1,m1+1,l2-1,m2+1,LLoad0,C0,S0,w3js,Itab) ...                 
                )*sqrt(2*l1+1)*sqrt(2*(l2-1)+1)...
                    *sqrt(2-(m1==0))*sqrt(2-(m2==0));
              
            ImXdX(:,index)=...
                -(b1_l1m1*a1_l2m2*...
                    calcGaunt(l1-1,m1-1,l2,m2-1,LLoad0,C0,S0,w3js,Itab)+...                     
                  b1_l1m1*a2_l2m2*...
                    calcGaunt(l1-1,m1-1,l2,m2+1,LLoad0,C0,S0,w3js,Itab)+...
                  b2_l1m1*a1_l2m2*...
                    calcGaunt(l1-1,m1+1,l2,m2-1,LLoad0,C0,S0,w3js,Itab)+...
                  b2_l1m1*a2_l2m2*...
                    calcGaunt(l1-1,m1+1,l2,m2+1,LLoad0,C0,S0,w3js,Itab) ...                 
                )*sqrt(2*(l1-1)+1)*sqrt(2*l2+1)...
                    *sqrt(2-(m1==0))*sqrt(2-(m2==0));
                            
                waitbar(2*index/((lenm^2)+lenm),h);
           end
       end
       delete(h)
                     
       % NOTE: The following is copied mutatis mutandis from the gauss
       % legendre case
       % In our ordering, the -1 precedes 1 and stands for the cosine term 
       comdex=[1:((lenm^2)+lenm)/2]';
       coss=gamini(comdex,comdi);
       % Need a vector of length "index" that points to the right
       % combination in XlmXlmp for the next array we are
       % designing. First, find the positions we've been missing 
       h=[dimK:-1:1']'; k=find(dems); kk=k+[1:length(k)]';
       % Where to insert other elements
       inpo=[indeks(cumsum(skip(h,kk)),k)+1]';
       % How many elements to insert
       inel=h(kk);
       % Which elements to insert
       beg=inpo-h(k)+[1:length(inel)]';
       ent=inpo-h(k)+inel+[0:length(inel)-1]';
       ins=[];
       for ind=1:length(beg)
           ins=[ins coss(beg(ind):ent(ind))];
       end           
       thvals=thvals(:);                       
       
       % And how to do it
       bigo=insert(coss,ins,gamini(inpo,inel));
       % Get the longitudinal integration info for the domain
       if isstr(dom)
	      switch dom
	      case 'patch'
               % Get the parameters of the dom
               phint=dphpatch(acos(thvals),thR,th0,ph0);
	      case 'sqpatch'
               % Always the same longitudinal integration interval
               phint=repmat([phW phE],length(thvals),1);

          otherwise
               defval('Nk',10);
	      % Now we may have multiple pairs
	      phint=dphregion(acos(thvals)*180/pi,Nk,dom); 
	      phint=phint*pi/180;
	      end
       else
	   % Now we may have multiple pairs
	   phint=dphregion(acos(thvals)*180/pi,[],dom); 
	   phint=phint*pi/180;
       end
 
       % Calculate the integral differences to get integration between 
       % two theta points           
       diffIdXdX=IdXdX(1:end-1,:)-IdXdX(2:end,:);
       diffImXmX=ImXmX(1:end-1,:)-ImXmX(2:end,:);
       diffIdXmX=IdXmX(1:end-1,:)-IdXmX(2:end,:);
       diffImXdX=ImXdX(1:end-1,:)-ImXdX(2:end,:);  
       
       % The number of elements that will be calculated is
       nel=(dimK^2+dimK)/2;
       % Calculate matrix elements
       index=0;
       undex=0;
       andex=1;
       h=waitbar(0,'KERNELB: Evaluating integrals and assembling matrix');
       for lm1dex=1:dimK
           l1=bigl(lm1dex);
           m1=bigm(lm1dex);
           index=index+1;
           ondex=0;
           I =repmat(NaN,nstripes+1,dimK-lm1dex+1);
          dI1=repmat(NaN,nstripes+1,dimK-lm1dex+1); 
          dI2=repmat(NaN,nstripes+1,dimK-lm1dex+1); 
         ddI =repmat(NaN,nstripes+1,dimK-lm1dex+1); 
           for lm2dex=lm1dex:dimK
               l2=bigl(lm2dex);
               m2=bigm(lm2dex);
               ondex=ondex+1;
               undex=undex+1;
               waitbar(undex/nel,h);             
               if m1>0 & m2>0
                    I (:,ondex)= sinsin(acos(thvals),m1,m2,phint);
                   dI1(:,ondex)= sincos(acos(thvals),m1,m2,phint);
                   dI2(:,ondex)= sincos(acos(thvals),m2,m1,phint);
                  ddI (:,ondex)= coscos(acos(thvals),m1,m2,phint);              
               elseif m1<=0 & m2<=0
                    I (:,ondex)= coscos(acos(thvals),m1,m2,phint);
                   dI1(:,ondex)= sincos(acos(thvals),m2,m1,phint);
                   dI2(:,ondex)= sincos(acos(thvals),m1,m2,phint);
                  ddI (:,ondex)= sinsin(acos(thvals),m1,m2,phint);
               elseif m1>0 & m2<=0 
                    I (:,ondex)= sincos(acos(thvals),m1,m2,phint);
                   dI1(:,ondex)= sinsin(acos(thvals),m1,m2,phint); 
                   dI2(:,ondex)= coscos(acos(thvals),m1,m2,phint);
                  ddI (:,ondex)= sincos(acos(thvals),m2,m1,phint); 
               elseif m1<=0 & m2>0
                    I (:,ondex)= sincos(acos(thvals),m2,m1,phint);
                   dI1(:,ondex)= coscos(acos(thvals),m1,m2,phint);
                   dI2(:,ondex)= sinsin(acos(thvals),m1,m2,phint);
                  ddI (:,ondex)= sincos(acos(thvals),m1,m2,phint); 
               end     

                  I (:,ondex)=  I (:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1));
                 dI1(:,ondex)= dI1(:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1)); 
                 dI2(:,ondex)= dI2(:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1));
                ddI (:,ondex)=ddI (:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1));  
                
           end
	       % Assemble the kernel.
           % We want to integrate \int_thS^thN XlmXlmp(th)* 
           % \int_{phiW(th)}^phiE(th) E(phi) dphi dth. The problem is that 
           % by the shape of the continent, phiW and phiE depend on th. we
           % are pulling the phi integral out of the theta integral by the
           % trapezoidal rule.
           
           % Use a trapezoidal rule approach to get the phi integral
           % independent of the theta integral           
           trapI  =0.5*(I(1:end-1,:)+I(2:end,:));    
           trapdI1=0.5*(dI1(1:end-1,:)+dI1(2:end,:));  
           trapdI2=0.5*(dI2(1:end-1,:)+dI2(2:end,:));  
           trapddI=0.5*(ddI(1:end-1,:)+ddI(2:end,:));  
  
          
           B(index,lm1dex:dimK)=...               
               sum(trapI  .*diffIdXdX(:,bigo(andex:undex)),1) + ...
               sum(trapddI.*diffImXmX(:,bigo(andex:undex)),1);
           
           D(index,lm1dex:dimK)=...
               sum(trapdI1.*diffIdXmX(:,bigo(andex:undex)),1) + ...
               sum(trapdI2.*diffImXdX(:,bigo(andex:undex)),1); 

           % And symmetrize them
           B(lm1dex+1:dimK,index)= B(index,lm1dex+1:dimK)';
           D(lm1dex+1:dimK,index)=-D(index,lm1dex+1:dimK)';

           andex=undex+1;
       end
        % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
        % did above here, taking the output of YLM and multiplying
        B=B/4/pi;
        D=D/4/pi;

        delete(h)
        % Save this now
        if isstr(pars)
           dom=pars;
        end                                              
    
    elseif ~isstr(method) % GL with 'method' Gauss points
        intv=cos([thS thN]);
        nGL=max(ngl,2*Lmax);

        [w,x,N]=gausslegendrecof(nGL,[],intv);
        disp(sprintf('%i Gauss-Legendre points and weights calculated',N))

        % First calculate the Legendre functions themselves
        % Note that dimK==sum(dubs)
        dubs=repmat(2,lenm,1); dubs(mz)=1; comdi=[];     
        Xlm=repmat(NaN,length(x),lenm);
        dXlm=repmat(NaN,length(x),lenm);

        % Calculate the Legendre polynomials        
        ind=0;
        for l=0:Lmax     
            [X dX]=libbrecht(l,x(:)','sch',[]);  
             Xlm(:,ind+1:ind+l+1)=( X*sqrt(2*l+1))';
            dXlm(:,ind+1:ind+l+1)=(dX*sqrt(2*l+1))';                                                
            ind=ind+l+1;
        end

        h=waitbar(0,'KERNELB: Calculating all Legendre products');
        product1=repmat(NaN,length(x),((lenm^2)+lenm)/2);
        product2=repmat(NaN,length(x),((lenm^2)+lenm)/2);
        product3=repmat(NaN,length(x),((lenm^2)+lenm)/2);
        product4=repmat(NaN,length(x),((lenm^2)+lenm)/2);

        index=0;
        for lm1dex=1:lenm
            l1=dels(lm1dex);
            m1=dems(lm1dex);
            pos1=1/2*(l1)^2+1/2*l1+m1+1;
            comdi=[comdi ; dubs(lm1dex:lenm)];
            % Note that the last index will be ((lenm^2)+lenm)/2
            for lm2dex=lm1dex:lenm
                l2=dels(lm2dex);
                m2=dems(lm2dex);
                index=index+1;
                pos2=1/2*(l2)^2+1/2*l2+m2+1;
                % Calculate products of Legendre functions
                product1(:,index)= dXlm(:,pos1).*dXlm(:,pos2); 
                product2(:,index)= (m1*Xlm(:,pos1)).*(m2*Xlm(:,pos2))./...
                    (sin(acos(x)).*sin(acos(x)));  
                product3(:,index)= dXlm(:,pos1).*(m2*Xlm(:,pos2))./...
                    sin(acos(x));
                product4(:,index)=-(m1*Xlm(:,pos1)./sin(acos(x))).*...
                    dXlm(:,pos2);
                waitbar(2*index/((lenm^2)+lenm),h);
            end
        end
        delete(h)        

        % Note: Xlmlmp is length ((lenm^2)+lenm)/2 because the Legendre
        % products have a redundant half (almost), and can be ordered 
        % as [0 11 222].  In order to use this with the kernel, which is
        % length (dimK^2+dimK)/2, we need an indexing array of the same
        % length which is filled with indices to Xlmlmp.  This array (bigo)
        % fills Xlmlmp back out to the ordering used in the kernel, [0 111
        % 22222].  Each Legendre polynomial has the shortened ordering, so
        % Xlmlmp essentially has redundancy in two dimensions.  When "coss"
        % is made, this will expand Xlmlmp in basically one dimension.  The
        % second pass is made when "ins" is inserted at certain points into
        % "coss" to form "bigo."  Later, the Legendre products will be 
        % multiplied by different "I" matrices, representing the sine and 
        % cosine products  from the longitudinal integrals. In order to use 
        % a similar implementation as in the radial case, the kernel 
        % entries for l=0 are also calculated (=NaN) and later removed
        % 
        % For more information on this or other functions, see the 
        % Simons' group wiki page.

        % In our ordering, the -1 precedes 1 and stands for the cosine term 
        comdex=[1:((lenm^2)+lenm)/2]';
        coss=gamini(comdex,comdi);
        % Need a vector of length "index" that points to the right
        % combination in XlmXlmp for the next array we are
        % designing. First, find the positions we've been missing 
        h=[dimK:-1:1']'; k=find(dems); kk=k+[1:length(k)]';
        % Where to insert other elements
        inpo=[indeks(cumsum(skip(h,kk)),k)+1]';
        % How many elements to insert
        inel=h(kk);
        % Which elements to insert
        beg=inpo-h(k)+[1:length(inel)]';
        ent=inpo-h(k)+inel+[0:length(inel)-1]';
        ins=[];
        for ind=1:length(beg)
            ins=[ins coss(beg(ind):ent(ind))];
        end    
        % And how to do it
        bigo=insert(coss,ins,gamini(inpo,inel));
        % Get the longitudinal integration info for the domain    
        if isstr(dom)
            switch dom
                case 'patch'
                    % Get the parameters of the dom
                    phint=dphpatch(acos(x),thR,th0,ph0);
                case 'sqpatch'
                    % Always the same longitudinal integration interval
                    phint=repmat([phW phE],length(x),1);
                otherwise
                    defval('Nk',10);
                    % Now we may have multiple pairs
                    phint=dphregion(acos(x)*180/pi,Nk,dom);
                    phint=phint*pi/180;
            end
        else
            % Now we may have multiple pairs
            phint=dphregion(acos(x)*180/pi,[],dom);
            phint=phint*pi/180;
        end 

        % The number of elements that will be calculated is
        nel=(dimK^2+dimK)/2;
        % Calculate matrix elements
        index=0;
        undex=0;
        andex=1;
        h=waitbar(0,'KERNELB: Evaluating integrals and assembling matrix');
        for lm1dex=1:dimK
            l1=bigl(lm1dex);
            m1=bigm(lm1dex);
            index=index+1;
            ondex=0;
              I =repmat(NaN,length(x),dimK-lm1dex+1);
             dI1=repmat(NaN,length(x),dimK-lm1dex+1);
             dI2=repmat(NaN,length(x),dimK-lm1dex+1);
            ddI =repmat(NaN,length(x),dimK-lm1dex+1);
            for lm2dex=lm1dex:dimK
                l2=bigl(lm2dex);
                m2=bigm(lm2dex);
                ondex=ondex+1;
                undex=undex+1;
                waitbar(undex/nel,h);
                % Now evaluate the longitudinal integrals at the GL points
                % Important: Here, the negative m inside sin/cos must be
                % taken into account! This affects all the sin with a
                % negative m
                if m1>0 & m2>0
                      I (:,ondex)= sinsin(acos(x),m1,m2,phint);
                     dI1(:,ondex)= sincos(acos(x),m1,m2,phint);
                     dI2(:,ondex)= sincos(acos(x),m2,m1,phint);
                    ddI (:,ondex)= coscos(acos(x),m1,m2,phint);
                elseif m1<=0 & m2<=0
                      I (:,ondex)= coscos(acos(x),m1,m2,phint);
                     dI1(:,ondex)= sincos(acos(x),m2,m1,phint); % -m
                     dI2(:,ondex)= sincos(acos(x),m1,m2,phint); % -m
                    ddI (:,ondex)= sinsin(acos(x),m1,m2,phint);
                elseif m1>0 & m2<=0 % Got rid of redundant ,pars below here
                      I (:,ondex)= sincos(acos(x),m1,m2,phint);
                     dI1(:,ondex)= sinsin(acos(x),m1,m2,phint); % -m
                     dI2(:,ondex)= coscos(acos(x),m1,m2,phint);
                    ddI (:,ondex)= sincos(acos(x),m2,m1,phint); % -m
                elseif m1<=0 & m2>0
                      I (:,ondex)= sincos(acos(x),m2,m1,phint);
                     dI1(:,ondex)= coscos(acos(x),m1,m2,phint);
                     dI2(:,ondex)= sinsin(acos(x),m1,m2,phint); % -m
                    ddI (:,ondex)= sincos(acos(x),m1,m2,phint); % -m
                end                   

                  I (:,ondex)=  I (:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1));
                 dI1(:,ondex)= dI1(:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1)); 
                 dI2(:,ondex)= dI2(:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1));
                ddI (:,ondex)=ddI (:,ondex)/sqrt(l1*(l1+1)*l2*(l2+1));                               
            end
       
            B(index,lm1dex:dimK)=(w(:)'*(...                                    
                                     product1(:,bigo(andex:undex)).*I  +...
                                     product2(:,bigo(andex:undex)).*ddI ...                                     
                                                  ));   

            D(index,lm1dex:dimK)=(w(:)'*(...
                                     product3(:,bigo(andex:undex)).*dI1+...
                                     product4(:,bigo(andex:undex)).*dI2 ...
                                                  ));  

            % And symmetrize them           
            B(lm1dex+1:dimK,index)= B(index,lm1dex+1:dimK)';  
            D(lm1dex+1:dimK,index)=-D(index,lm1dex+1:dimK)';  
            andex=undex+1;                                                                                               
        end 
        delete(h)

        % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
        % did above here, taking the output of YLM and multiplying
        B=B/4/pi;
        D=D/4/pi;
        
    end
              
    % By whichever way you obtained the kernel, now check if you might want
    % to rotate it back so its eigenvectors are "right", right away, for
    % Antarctica without needing to rotate as part of LOCALIZATION
    if rotb==1
      disp('The input coordinates were rotated. Kernel will be unrotated,')
      disp('so its eigenfunctions will show up in the right place')
      disp(' ')
      % Get the rotation parameters for this particular region
      [XY,lonc,latc]=eval(sprintf('%s(%i)',dom,pars));
      % Getting rid of the NaNs in the extended B and D or there
      
      B(:,1)=zeros(size(B,1),1);
      B(1,:)=zeros(1,size(B,2));
      D(:,1)=zeros(size(D,1),1);
      D(1,:)=zeros(1,size(D,2));
      
      % Rotate the kernels. Why the transposed for D?  
      B=klmlmp2rot(B,lonc,latc);
      D=klmlmp2rot(D',lonc,latc);      
          
    else
      [lonc,latc]=deal(0);
    end
   
    % because the vector Slepian horizontal components only start with
    % degree l=1, remove the first row and column
    B=B(2:end,2:end);
    D=D(2:end,2:end);
    K=[B D;D' B];
    
    if rotb==1
      disp(sprintf('Asymmetry of kernel through rotation is %g',...
		   norm(K-K')));
      disp('Getting rid of this by storing (K+K^T)/2')
      K=(K+K')/2;
    end
    
    % Save this now
    if isstr(pars)
        dom=pars;
    end
    
    save(fnpl,'Lmax','K','B','D','dom','ngl','XY',...
   'lonc','latc')

    end % If loading or calculating
end

varns={K,B,D,XY};
varargout=varns(1:nargout);
 
elseif strcmp(Lmax,'demo1')
    dom='africa'%'australia'
    Lmax=12;
    ngl=500;  
    P=kernelc(Lmax,dom,[],ngl); 
    K=kernelb(Lmax,dom,[],ngl);
    area=P(1,1)*4*pi;
    shannon=(length(P)+length(K))*area/(4*pi);
    closeup=ceil(3*shannon);
    EP=sort(eig(P),'descend');
    EK=sort(eig(K),'descend');
    E=sort([EP;EK],'descend');
    %plot(E(1:closeup),'-xk')
    plot(E,'-xk')
    title('Eigenvalues')
    hold on    
    %plot(EP(1:closeup),':r')
    %plot(EK(1:closeup),':')
    plot(EP,':r')
    plot(EK,':')
    plot([shannon shannon],[0 1],'--k')
    legend('Total','Radial','Tangential','Shannon')
    hold off    
    figdisp('Eigenspec',sprintf('%s_%i',dom,Lmax))  

elseif strcmp(Lmax,'demo2')
    dom='africa';
    Lmax=5;
    tic
    system(sprintf('rm ./KERNELB/WREG-%s-%d-tang.mat',dom,Lmax))
    [Kgl,Bgl,Dgl]=kernelb(Lmax,dom,[]);
    timegl=toc;
    disp(sprintf('Time for GL: %g sec',timegl))
    tic
    system(sprintf('rm ./KERNELB/WREG-%s-%d-tang.mat',dom,Lmax))
    [Kpaul,Bpaul,Dpaul]=kernelb(Lmax,dom,[],'paul',[]);
    timepaul=toc;
    disp(sprintf('Time for paul: %g sec',timepaul))     
    maxB=max(max(abs((Bpaul-Bgl))));
    maxD=max(max(abs((Dpaul-Dgl))));
    disp(sprintf('Max difference in matrix entries %d',...
        max(maxB,maxD)));
   
    Epaul=sort(eig(Kpaul),'descend');
    Egl=sort(eig(Kgl),'descend');
    
    plot(Epaul,'--r')
    hold on
    plot(Egl)
    hold off
    legend('Eigenvalues Paul','Eigenvalues Gauss-Legendre')


    
end




function varargout=kernelbp(Lmax,dom,pars,method,rotb)
% [K,B,D,XY]=KERNELBP(Lmax,dom,pars,method,rotb)
%  
% Parallel version of the vector spherical harmonic localization kernel for
% the tangential plane for a all degrees l between 0 and Lmax and for all
% orders -m<=l<=m. Returns the full tangential Kernel and the B = C and D 
% component.
%
% INPUT:
%
% Lmax      Bandwidth
% dom       Either the region 'england', 'eurasia', 'australia',
%           'greenland', 'africa', 'samerica', 'amazon', 'orinoco', 
%           'gpsnamerica', 'antarctica', 'alloceans', 'namerica' [default]
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
%           [default: 10]
%           OR: the string with the name you want the result saved as
% method    an interger for Gauss-Legendre with the indicated number of
%           Gauss-Legendre points or 'paul' for Paul-Gaunt
% rotb      0 That's it, you got it [default: 0]
%           1 For, e.g., 'antarctica', if you were given rotated 
%           coordinates to make the integration procedure work, this option 
%           makessure that the kernel matrix reflects this. If not, you 
%           have to apply counterrotation after diagonalizing in
%           LOCALIZATION.
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
% kernelbp('demo1') Plot the eigenvalues for Australia
%
% kernelbp('demo2') Compare the performance of the parallel to the serial
%                   version of kernelb
%
% kernelbp('demo3') Compare Antarctica "first rotating the kernel and then
%                   finding the eigenfields" to "first finding the
%                   eigenfields and then rotate them"
%
% kernelbp('demo4') Calculate eigenvalues of a "belt"-region using the
%                   analytic calculation KERNALTANCAPM and the numerical
%                   calculation KERNELBP and compare them
%
% kernelbp('demo5') Plot and compare Slepians for a "belt"-region using the
%                   analytic and the numerical kernel
%
% kernelbp('demo6') Calculate eigenvalues of a spherical cap using the
%                   analytic kernel with rotation and using the numerical 
%                   kernel  
%
% kernelbp('demo7') plot and compare Slepians for a spherical cap region
%                   using the analytic and the numerical kernel 
%
% On 06/29/2017, plattner-at-alumni.ethz.ch made
% Antartica rotate back by default
%  
% Last modified by plattner-at-alumni.ethz.ch, 07/14/2017
% 
% See also VECTORSLEPIAN, BLMCLM2XYZ, KERNELB, KERNELCP, KERNALTANCAPM
  
defval('Lmax',18)
defval('dom','namerica')
defval('method',200)
defval('rotb',0)

if strcmp(dom,'antarctica')
  rotb = 1;
end


if ~isstr(Lmax)
% Generic path name that I like
%filoc=fullfile(getenv('IFILES'),'KERNELBP');    
filoc=fullfile(getenv('IFILES'),'KERNELB');    


if ~(isstr(dom) || length(dom)>2)
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
    
    
    if exist(fnpl,'file')==2 && (~isstr(ngl))
        load(fnpl)
        disp(sprintf('%s loaded by KERNELBP',fnpl))
        K=[B D;D' B];      
        if rotb==1
            disp(sprintf('Asymmetry of kernel through rotation is %g',...
                norm(K-K')));
            disp('Getting rid of this by storing (K+K^T)/2')
            K=(K+K')/2;
        end
%         % Save this now
%         if isstr(pars)
%             dom=pars;
%         end
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
    if ~isstr(method) % GL with 'method' Gauss points       
        % See if we can run this calculation in parallel, and set a flag
        try
            parpool      
        end
        
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
             Xlm(:,ind+1:ind+l+1)=(X*sqrt(2*l+1))';
            dXlm(:,ind+1:ind+l+1)=(dX*sqrt(2*l+1))';                                                
            ind=ind+l+1;
        end

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
        % cosine products from the longitudinal integrals. In order to use
        % a similar implementation as in the radial case, the kernel 
        % entries for l=0 are also calculated (=NaN) and later removed
        % 
        % For more information on this or other functions, see the 
        % Simons' group wiki page.

        % In our ordering, the -1 precedes 1 and stands for the cosine term         
        % comdex=[1:((lenm^2)+lenm)/2]';
        % coss=gamini(comdex,comdi);
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
        parfor lm1dex=1:dimK
            l1=bigl(lm1dex);
            m1=bigm(lm1dex);
            % Can only use the loop variable once per index. So for Klmlmp,
            % also use index=lm1dex
            index=lm1dex;          
            ondex=0;
              I =repmat(NaN,length(x),dimK-index+1);
             dI1=repmat(NaN,length(x),dimK-index+1);
             dI2=repmat(NaN,length(x),dimK-index+1);
            ddI =repmat(NaN,length(x),dimK-index+1);
            % Instead of counting up andex and undex, write expressions for
            % them analytically for the row of B and D we are calculating
            countdown = [dimK:-1:1];
            andex = 1 + sum(countdown(1:(index-1)));
            undex = sum(countdown(1:index));
            % We know bigo, andex, and undex, so just calculated exactly 
            % which parts of XlmXlmp you need for this specific iteration
            smalll1 = abs(l1);
            smallm1 = abs(m1);
            pos1=1/2*(smalll1)^2+1/2*smalll1+smallm1+1;
            product1=0; product2=0; product3=0; product4=0;
            for lm2dex=lm1dex:dimK
                l2=bigl(lm2dex);
                m2=bigm(lm2dex);
                smalll2 = abs(l2);
                smallm2 = abs(m2);
                pos2=1/2*(smalll2)^2+1/2*smalll2+smallm2+1;
                ondex=ondex+1;
                product1(1:length(x),ondex)= dXlm(:,pos1).*dXlm(:,pos2);
                product2(1:length(x),ondex)= (smallm1*Xlm(:,pos1)).*...
                    (smallm2*Xlm(:,pos2))./(sin(acos(x)).*sin(acos(x)));
                product3(1:length(x),ondex)= dXlm(:,pos1).*...
                    (smallm2*Xlm(:,pos2))./sin(acos(x));
                product4(1:length(x),ondex)=-(smallm1*Xlm(:,pos1)./...
                    sin(acos(x))).*dXlm(:,pos2);                                               
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

            % Do the calculation and set as a temp variable                      
                        
            temprowB=(w(:)'*(...                                    
                             product1.*I + ...
                             product2.*ddI ...                                     
                            ));   
                        
            temprowD=(w(:)'*(...
                             product3.*dI1 +...
                             product4.*dI2  ...
                            ));                                                                         
                        
            % Pad the temp variable with the appropriate zeros out front
            temprowB = [zeros(1,(index-1)) temprowB];
            temprowD = [zeros(1,(index-1)) temprowD];          
             
            % Now we can distribute over the kernel.  We need to do it this 
            % way because if you slice Klmlmp with the loop variable 
            % (lm1dex) then all other indicies need to be constant, 
            % or ':', or 'end.'
            B(lm1dex,:)=temprowB;
            D(lm1dex,:)=temprowD;
        end %parfor
            
        % Close the parpool
        delete(gcp('nocreate'))
      
        % Symmetrize the Kernel
        B = B + B' - diag(diag(B));
        D = D - D';

        % To make this exactly equivalent to Tony's \ylm, i.e. undo what we
        % did above here, taking the output of YLM and multiplying
        B=B/4/pi;
        D=D/4/pi;
             
        % By whichever way you obtained the kernel, now check if you might 
        % want to rotate it back so its eigenvectors are "right", right 
        % away, for Antarctica without needing to rotate as part of 
        % LOCALIZATION
    if rotb==1
      disp('The input coordinates were rotated. Kernel will be unrotated,')
      disp('so its eigenfunctions will show up in the right place')
      disp(' ')
      % Get the rotation parameters for this particular region
      [XY,lonc,latc]=eval(sprintf('%s(%i)',dom,pars));
      % Get rid of the NaNs in the extended B and D or there
      B(:,1)=zeros(size(B,1),1);
      B(1,:)=zeros(1,size(B,2));
      D(:,1)=zeros(size(D,1),1);
      D(1,:)=zeros(1,size(D,2));
            
      % Rotate the kernels, properly. Why the transposed for D?                                       
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
        disp('Get rid of this by storing (K+K^T)/2')
        K=(K+K')/2;
    end
    % Save this now
    if isstr(pars)
        dom=pars;
    end


    save(fnpl,'Lmax','B','D','dom','ngl','XY',...
         'lonc','latc')
 end      
    end % If loading or calculating
end

varns={K,B,D,XY};
varargout=varns(1:nargout);

 
elseif strcmp(Lmax,'demo1')
    dom='australia'
    Lmax=22;
    ngl=500; 
    P=kernelcp(Lmax,dom,[],ngl);  
    K=kernelbp(Lmax,dom,[],ngl);
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
    % Comparing the calculation speed between kernelb and Kernelbp
    dom='australia'
    Lmax=10;
    ngl=500; 
    system(sprintf('rm ./KERNELBP/WREG-australia-%d-tang.mat',Lmax))
    system(sprintf('rm ./KERNELB/WREG-australia-%d-tang.mat',Lmax))
    tic;
    [K,B,D]=kernelb(Lmax,dom,[],ngl);
    time=toc;
    disp(sprintf('Time for serial: %g sec',time));
    tic;
    [Kp,Bp,Dp]=kernelbp(Lmax,dom,[],ngl);
    timep=toc;
    disp(sprintf('Time for parallel: %g sec',timep));    
    disp(sprintf('norm of serial: %g, norm of difference: %g',...
        norm(K),norm(K-Kp)));
    
elseif strcmp(Lmax,'demo3')
    % Plotting a tangential vector Slepian for Antarctica
    % One time rotate the equatorial antarctica, other time calculate
    % rotated kernel
    clf;
    fig2print(gcf,'flandscape')
    dom='antarctica';
    comp='tangential';
    Lmax=18;
    index=5%10;
    res=[0.2 5];%[1 5]%[0.2 3];
    range=[0 360 -90 90-sqrt(eps)];  
    c11cmn=[range(1) range(4) range(2) range(3)];
    [~,lonc,latc]=eval(sprintf('%s(10)',dom));   
    [ah,ha,H]=krijetem(subnum(2,2));
    C=[]; V=[];
    
    % First rotate after calculating   
    rotb=0;
    [~,~,~,C1,V1,blmcosi1,clmcosi1]=vectorslepian(Lmax,dom,...
        comp,index,res,c11cmn,C,V,rotb);  
    % Rotate
    alp=-lonc;
    bta=latc;
    gam=0;
    % The rotation routine is written for plm coefficients. We must
    % therefore add the l=0 coeficients
    blmcosi1=[0 0 0 0;blmcosi1];
    clmcosi1=[0 0 0 0;clmcosi1];    
    % Because orthogonal operations commute with the operators that define
    % the vector spherical harmonics from the spherical harmonics, we can
    % rotate the coefficients using the rotation for spherical harmonics
    blmcosip=plm2rot(blmcosi1,alp,bta,gam); 
    % And delete the l=0 coefficient again
    blmcosip=blmcosip(2:end,:);  
    clmcosip=plm2rot(clmcosi1,alp,bta,gam);   
    % And delete the l=0 coefficient again also for the clm
    clmcosip=clmcosip(2:end,:);
    % Now calculate the data
    [datar{1},lon{1},lat{1}]=blmclm2xyz(blmcosip,clmcosip,res(1));
    % and on a less dense grid to show the vectors
    absdatar=sqrt(datar{1}(:,:,1).^2+datar{1}(:,:,2).^2);
    dmax=max(max(absdatar));
    [datar{2},lon{2},lat{2}]=blmclm2xyz(blmcosip,clmcosip,res(2));     
    axes(ah(1))
    % Plot the absolute values
    imagefnan([range(1) range(4)],[range(2) range(3)],absdatar,...
        kelicol,[-dmax dmax],[],1,100);
    hold on
    % And the directions
    quiverimage(datar{2},lon{2},lat{2})
    axis off
    plotcont;
    text(180,100,'Rotate after solving','HorizontalAlignment','center')
    hold off
    % Now plot the same on a three dimensional sphere
    axes(ah(3))
    % Plot the absolute value
    plotonearth(-absdatar)
    kelicol(1)
    caxis([-dmax dmax])    
    hold on
    % Plot a circle around the sphere such that the boundary is visible 
    circ;
    % Plot the directions
    quiversphere(datar{2},[],[],[],0.01)
    % Rotate the sphere to Antarctica
    view(90,-90)
    hold off
    axis off

    % Now rotate the kernel and then calculate the Slepian
    axes(ah(2))    
    rotb=1;
    % Because the same eigenvalue shows up twice (see paper), it is a
    % coincidence, which of the two mutually pointwise perpendicular
    % Slepian functions with the same eigenvalue shows up first. Hence it
    % might be necessary to compare index before with index+1 or index-1
    % here.
    %index=index-1;
    [data,lat,lon]=vectorslepian(Lmax,dom,comp,index,res,...
        c11cmn,C,V,rotb);  
    % For some reason the sign is different. This does not matter.
    data{2}=-data{2};
    absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
    dmax=max(max(absdata));
    imagefnan([range(1) range(4)],[range(2) range(3)],absdata,...
        kelicol,[-dmax dmax],[],1,100);
    hold on
    % Now the directions
    quiverimage(data{2},lon{2},lat{2})
    axis off
    % Plot the continents too
    plotcont;
    text(180,100,'Rotate before solving','HorizontalAlignment','center')
    hold off
    % Plot the same on a three dimensional sphere      
    axes(ah(4))
    plotonearth(-absdata)
    kelicol(1)
    caxis([-dmax dmax]) 
    hold on
    % Plot a circle around the sphere such that the boundary is visible 
    circ;
    % Plot the directions
    quiversphere(data{2},[],[],[],0.01)
    % Rotate the sphere to Antarctica
    view(90,-90)
    hold off
    axis off
     
    % cosmetics
    serre(ah(1:2),0.75,'across')
    serre(ah(3:4),0.75,'across')
    serre(ha(1:2),1.25,'down')
    serre(ha(3:4),1.25,'down')
    
    nrmdiff=min(sqrt(sum(sum(sum((data{1}-datar{1}).^2)))) , ...
                sqrt(sum(sum(sum((data{1}+datar{1}).^2)))) );
    nrm=sqrt(sum(sum(sum((datar{1}).^2))));
    
    disp(sprintf('Relaive difference vectors = %g',...
       nrmdiff/nrm));
   
    disp(sprintf('Relaive difference abs values = %g',...
       norm(absdatar-absdata)/norm(absdatar)));

    figdisp(comp,sprintf('%s_%i_%i_comparison',dom,Lmax,index))
    disp('Use the -r600 option')
    
elseif strcmp(Lmax,'demo4')
    % Compare squarepatch to polar cap
    THN=30;
    THS=180-THN;
    Lmax=10;
    thS=THS*pi/180;
    thN=THN*pi/180;
    pars=[thN thS -pi pi];
    K1=kernelbp(Lmax,'sqpatch',pars);    
    
    m=0;
    K=kerneltancapm([THS THN],Lmax,m); 
    for m=1:Lmax
        [Mm,Mmm]=kerneltancapm([THS THN],Lmax,m);
        K=blkdiag(K,Mm,Mmm);
    end      
    Ebm=eig(K);
    Eb=sort(Ebm,'descend');
    plot(Eb,'x')
    hold on
    Esqpatch=eig(K1);
    Es=sort(Esqpatch,'descend');
    plot(Es,'k--')
    hold off
    
elseif strcmp(Lmax,'demo5')
    % Plot sqpatch Slepian 
    % Compare squarepatch to polar cap
    THN=80;
    THS=110;%180-THN;
    thS=THS*pi/180;
    thN=THN*pi/180;
    Lmax=10;
    index=1;
    res=[1 8];
    pars=[thN thS 0 2*pi];%[thN thS 0 2*pi-eps];
    range=[0 360 -90 90-sqrt(eps)];  
    c11cmn=[range(1) range(4) range(2) range(3)];
    
    K=kernelbp(Lmax,'sqpatch',pars);
    [C,V]=eig(K);    
    [V,isrt]=sort(sum(V,1),'descend');
    C=C(:,isrt(1:length(K)));    
    [blmcosi,clmcosi]=coef2blmclm(C(:,index),Lmax);
    
    [ah,ha,H]=krijetem(subnum(1,2));
    axes(ah(1))
    [datan{1},lon,lat]=blmclm2xyz(blmcosi,clmcosi,res(1),c11cmn);
    [datan{2},lon,lat]=blmclm2xyz(blmcosi,clmcosi,res(2),c11cmn);
    absdata1=sqrt(datan{1}(:,:,1).^2+datan{1}(:,:,2).^2);
    dmax=max(max(absdata1));
    imagefnan([range(1) range(4)],[range(2) range(3)],absdata1,...
        kelicol,[-dmax dmax],[],1,100);
    title(sprintf('%s=%1.9f','\lambda',V(index)));
    hold on
    quiverimage(datan{2});
    hold off
    
    [data,lat,lon,C,V,Vtot]=capvectorslepian(Lmax,[THS THN],[],...
        index,res,[],[],[],c11cmn); 
    axes(ah(2))
    absdata2=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
    dmax=max(max(absdata2));
    imagefnan([range(1) range(4)],[range(2) range(3)],absdata2,...
        kelicol,[-dmax dmax],[],1,100);
    title(sprintf('%s_{%i} =%1.9f; m = %i',...
           '\lambda',Vtot(index,3),Vtot(index,1),Vtot(index,2)));
    hold on
    quiverimage(data{2});
    hold off    
    nrmdiff=min(sqrt(sum(sum(sum((datan{1}-data{1}).^2)))) , ...
                sqrt(sum(sum(sum((datan{1}+data{1}).^2)))) );
    nrm=sqrt(sum(sum(sum((data{1}).^2))));
    
    disp(sprintf('Relaive difference vectors = %g',...
       nrmdiff/nrm));
   
    disp(sprintf('Relaive difference abs values = %g',...
       norm(absdata1-absdata2)/norm(absdata1)));
   
    
elseif strcmp(Lmax,'demo6')    
    TH=20;
    th0=pi/2;
    ph0=0;
    Lmax=15;
    pars=[th0 ph0 TH/180*pi];
    K=kernelbp(Lmax,'patch',pars);
    [Cnum,Vnum]=eig(K);    
    [Vnum,isrt]=sort(sum(Vnum,1),'descend');
    
    [data,lat,lon,C,V,Vtot]=capvectorslepian(Lmax,TH);
    
    plot(Vnum,'rx');
    hold on
    plot(Vtot(:,1),'o')
    
    rms=sqrt(sum((Vtot(:,1)-Vnum').^2)/length(Vnum));
    disp(sprintf('Eigenvalue rms difference is %g',rms));

elseif strcmp(Lmax,'demo7')
    % Plot the index best Slepian for rotated and for numerically
    % calculated cap over th0 ph0 with opening angle TH
    TH=40;
    th0=pi/2;
    ph0=pi;
    Lmax=18; 
    index=5;    
    % Plotting stuff
    res=[1 5];
    range=[0 360 -90 90-sqrt(eps)];  
    c11cmn=[range(1) range(4) range(2) range(3)];
    [ah,ha,H]=krijetem(subnum(2,1));
    
    % Fist the numerical solution
    axes(ah(1))
    pars=[th0 ph0 TH/180*pi];
    % Solve the Slepian problem for this spherical cap numerically
    K=kernelbp(Lmax,'patch',pars);
    [Cnum,Vnum]=eig(K);    
    [Vnum,isrt]=sort(sum(Vnum,1),'descend');
    Cnum=Cnum(:,isrt(1:length(K)));    
    % Turn the Slepian coefficients into pointwise data
    [blmcosi,clmcosi]=coef2blmclm(Cnum(:,index),Lmax);   
    [datan{1},lon,lat]=blmclm2xyz(blmcosi,clmcosi,res(1),c11cmn);
    [datan{2},lon,lat]=blmclm2xyz(blmcosi,clmcosi,res(2),c11cmn);
    % Plot the pointwise data
    absdatan1=sqrt(datan{1}(:,:,1).^2+datan{1}(:,:,2).^2);
    dmax=max(max(absdatan1));
    imagefnan([range(1) range(4)],[range(2) range(3)],...
       absdatan1,kelicol,[-dmax dmax],[],1,100);
    title(sprintf('Numerical, %s=%1.9f','\lambda',Vnum(index)));
    hold on
    quiverimage(datan{2});
    circ(TH,[],[ph0 pi/2-th0]*180/pi,[],'LineStyle','--')
    hold off
    
    % Now the analytical solution including rotation    
    axes(ah(2))
    % Solve the Slepian problem for the polar spherical cap analyically
    [~,~,~,C,V,Vtot,blmcosia,clmcosia]=capvectorslepian(Lmax,TH,[],index);
    % Now rotate to the location chosen for the numerical solution
    alp=0;
    bta=th0*180/pi;
    gam=ph0*180/pi-180;
    blmcosip=plm2rot([0 0 0 0;blmcosia],alp,bta,gam);
    clmcosip=plm2rot([0 0 0 0;clmcosia],alp,bta,gam);
    blmcosip=blmcosip(2:end,:);
    clmcosip=clmcosip(2:end,:);
    % Turn the Slepian coefficients into pointwise data
    [dataa{1},lon,lat]=blmclm2xyz(blmcosip,clmcosip,...
        res(1),c11cmn);
    [dataa{2},lon,lat]=blmclm2xyz(blmcosip,clmcosip,...
        res(2),c11cmn);
    % Plot the pointwise data
    absdataa1=sqrt(dataa{1}(:,:,1).^2+dataa{1}(:,:,2).^2);
    dmax=max(max(absdataa1));
    imagefnan([range(1) range(4)],[range(2) range(3)],...
       absdataa1,kelicol,[-dmax dmax],[],1,100);
    title(sprintf('Analytical, %s=%1.9f','\lambda',Vtot(index,1)));
    hold on
    quiverimage(dataa{2});
    circ(TH,[],[ph0 pi/2-th0]*180/pi,[],'LineStyle','--')
    hold off

    % Display error
    rmsdiff1=sqrt(sum(sum(sum((dataa{1}-datan{1}).^2)))...
        /size(dataa{1},1)/size(dataa{1},2)/size(dataa{1},3));
    rmsdiff2=sqrt(sum(sum(sum((dataa{1}+datan{1}).^2)))...
        /size(dataa{1},1)/size(dataa{1},2)/size(dataa{1},3));
    rms=sqrt(sum(sum(sum((dataa{1}).^2)))...
        /size(dataa{1},1)/size(dataa{1},2)/size(dataa{1},3));  
    disp(sprintf(...
        'Rms difference = %g, rms value of analytical field = %g'...
        ,min(rmsdiff1, rmsdiff2),rms));
    
    % Display the coefficients
    figure
    [ah,ha,H]=krijetem(subnum(2,2));
    bcoeffmatrixn=intomatrix(blmcosi,Lmax,0);
    ccoeffmatrixn=intomatrix(clmcosi,Lmax,0);
    axes(ah(1))
    imagefnan([-Lmax,0],[Lmax,Lmax],bcoeffmatrixn);
    title('Numerical Blm coefficients')
    axes(ah(2))
    imagefnan([-Lmax,0],[Lmax,Lmax],ccoeffmatrixn);
    title('Numerical Clm coefficients')
    
    bcoeffmatrixa=intomatrix(blmcosip,Lmax,0);
    ccoeffmatrixa=intomatrix(clmcosip,Lmax,0);
    axes(ah(3))
    % The coefficients can be multiplied with a -1 without changing the
    % optimality and normalization of the Slepian. Wich one is displayed is
    % random, therefore to have the same coefficient values sign, maybe
    % plot -bcoeffmatrixa and -ccoeffmatrixa instead of without the -
    imagefnan([-Lmax,0],[Lmax,Lmax],-bcoeffmatrixa);
    title('Analytic Blm coefficients')
    axes(ah(4))
    imagefnan([-Lmax,0],[Lmax,Lmax],-ccoeffmatrixa);
    title('Analytic Clm coefficients')
    
    disp('If they look different, there are two possibilities:')
    disp('1) The other eigenfunction with the same eigenvalue is displayed (see the plots of the fields)')
    disp('2) The coefficients have a different sign. Plot -coeffmatrices instead')       
    
end


function varargout=blmclm2xyz(blmcosi,clmcosi,degres,c11cmn,lmax)
% [r,lon,lat]=BLMCLM2XYZ(blmcosi,clmcosi,degres,c11cmn,lmax)
%
% Inverse tangential vector spherical harmonic transform
%
% Compute a tangental vector field (coefficents for theta and phi) from
% coefficients of the Blm and Clm given as [l m Ccos Csin] (not necessarily 
% starting from zero, but sorted) with degree resolution 'degres' 
% [default: approximate Nyquist degree].
%
% Using 4*pi-normalized real spherical harmonics
%
% INPUT:
%
% blmcosi   Matrix listing l,m,cosine and sine expansion coefficients for
%           Blm
% clmcosi   Matrix listing l,m,cosine and sine expansion coefficients for
%           Clm
% degres    Longitude/ latitude spacing, in degrees [default: Nyquist] OR
%           "lat": a column vector with latitudes [degrees]
% c11cmn    Corner nodes of lon/lat grid [default: 0 90-sqrt(eps) 360 -90]
%           OR "lon": a column vector with longitudes [degrees]
% lmax      Maximum bandwidth expanded at a time [default: 720]
%
% OUTPUT: 
%
% r         The field (matrix for a grid, vector for scattered points)
%           r(:,:,1) is the phi-component; r(:,:,2) is the theta-component.
%           First dimension of r is lat, second dimension is lon.
% lon,lat   The grid (matrix) or evaluation points (vector), in degrees
%
% See also KERNELB, KERNELBM
%
% Last modified by plattner-at-alumni.ethz.ch 01/25/2012

if ~isstr(blmcosi)
    % First check that the expansion of both, Blm and Clm are within the
    % same degrees
    if ((blmcosi(1)~=clmcosi(1))||(blmcosi(end,1)~=clmcosi(end,1)))
        error('Blm and Clm expansion are not for the same degrees!')
    end
            
    % Lowest degree of the expansion
    lmin=blmcosi(1);
    % Highest degree (bandwidth of the expansion)
    L=blmcosi(end,1);
    % Default resolution is the Nyquist degree; return equal sampling in
    % longitude and latitude; sqrt(L*(L+1)) is equivalent wavelength
    degN=180/sqrt(L*(L+1));
    defval('degres',degN);
    
    % When do you get a task bar?
    taskmax=100;
    
    % Default grid is all of the planet
    defval('c11cmn',[0 90-sqrt(eps) 360 -90]);
    
    % Build in maxima to save memory
    defval('latmax',Inf); 
    defval('lmax',999);
        
    % But is it a grid or are they merely scattered points?
    if length(degres)==length(c11cmn)
        % It's a bunch of points!
        nlat=length(degres);
        nlon=length(c11cmn);
        % Colatitude vector in radians
        theta=[90-degres(:)']*pi/180;
        % Longitude vector in radians
        phi=c11cmn(:)'*pi/180;

        % Initialize output vector
        r=zeros(nlat,1,2);
    
        % Now if this is too large reduce lmax, our only recourse to 
        % hardcode
        ntb=256;
        if round(sqrt(nlat)) >= ntb || round(sqrt(nlon)) >= ntb
            lmax=round(ntb);
        end
    elseif length(degres)==1 && length(c11cmn)==4
        % It's a grid
        if degres>degN
        disp('PLM2XYZ: You can do better! Ask for more spatial resolution')
        disp(sprintf('Spatial sampling ALLOWED: %8.3f ; REQUESTED: %6.3f',...
                degN,degres))
        end
        % The number of longitude and latitude grid points that will be 
        % computed
        nlon=min(ceil([c11cmn(3)-c11cmn(1)]/degres+1),latmax);
        nlat=min(ceil([c11cmn(2)-c11cmn(4)]/degres+1),2*latmax+1);

        % Initialize output grid
        r=zeros(nlat,nlon,2);
    
        % Longitude grid vector in radians
        phi=linspace(c11cmn(1)*pi/180,c11cmn(3)*pi/180,nlon);
        % Colatitude grid vector in radians
        theta=linspace([90-c11cmn(2)]*pi/180,[90-c11cmn(4)]*pi/180,nlat);
    else
        error('Make up your mind - is it a grid or a list of points?')
    end
    
    
    % Piecemeal degree ranges
    % Divide the degree range increments spaced such that the additional
    % number of degrees does not exceed addmup(lmax)
    % If this takes a long time, abort it
    els=0; ind=0;
    while els<L
        ind=ind+1;
        % Take positive root
        els(ind+1)=min(floor(max(roots(...
        [1 3 -els(ind)^2-3*els(ind)-2*addmup(lmax)]))),L);
        if any(diff(els)==0)
            error('Increase lmax as you are not making progress')
        end
    end
    % Now els contains the breakpoints of the degrees
    if ~all(diff(addmup(els))<=addmup(lmax))
        error('The subdivision of the degree scale went awry')
    end
    
    % Here's the lspacings
    if length(els)>2
        els=pauli(els,2)+...
        [0 0 ; ones(length(els)-2,1) zeros(length(els)-2,1)];
    end

    for ldeg=1:size(els,1)
        ldown=els(ldeg,1);
        lup=els(ldeg,2);
        % Evaluate Legendre polynomials at selected points
        try
            Plm=nan(length(theta),addmup(lup)-addmup(ldown-1));
        catch ME
            error(sprintf('\n %s \n\n Decrease lmax in PLM2XYZ \n',...
                ME.message))
        end
        if [lup-ldown]>taskmax && length(degres) > 9999
            h=waitbar(0,sprintf(...
                'Evaluating Legendre polynomials between %i and %i',...
                ldown,lup));
        end
        in1=0;
        in2=ldown+1;
        % Always start from the beginning in this array, regardless of lmin
        for l=ldown:lup
           [Plm dPlm]=libbrecht(l,cos(theta(:)'),'sch');
           divsin=repmat(1./sin(theta(:)'),size(Plm,1),1);        
           dphPlm(:,in1+1:in2)=(1/sqrt(l*(l+1))* Plm.*divsin*sqrt(2*l+1))';
           dthPlm(:,in1+1:in2)=(1/sqrt(l*(l+1))*dPlm        *sqrt(2*l+1))';
           in1=in2;
           in2=in1+l+2;
           if [lup-ldown]>taskmax && length(degres) > 9999
               waitbar((l-ldown+1)/(lup-ldown+1),h)
           end
        end
        if [lup-ldown]>taskmax && length(degres) > 9999
            delete(h)
        end
        
        % Loop over the degrees
        more off
        for l=max(lmin,ldown):lup
            % Compute Schmidt-normalized Legendre functions at 
            % the cosine of the colatitude (=sin(lat)) and 
            % renormalize them to the area of the unit sphere
      
            % Remember the Plm vector always starts from ldown
            b=addmup(l-1)+1-addmup(ldown-1);
            e=addmup(l)-addmup(ldown-1);

            dphplm=dphPlm(:,b:e)';
            dthplm=dthPlm(:,b:e)';
            
            m=0:l;
            mphi=m(:)*phi(:)';
      
    % Normalization of the harmonics is to 4\ pi, the area of the unit
    % sphere: $\int_{0}^{\pi}\int_{0}^{2\pi}
    % |P_l^m(\cos\theta)\cos(m\phi)|^2\sin\theta\,d\theta d\phi=4\pi$.
    % Note the |cos(m\phi)|^2 d\phi contributes exactly \pi for m~=0
    % and 2\pi for m==0 which explains the absence of sqrt(2) there;
    % that fully normalized Legendre polynomials integrate to 1/2/pi
    % that regular Legendre polynomials integrate to 2/(2l+1),
    % Schmidt polynomials to 4/(2l+1) for m>0 and 2/(2l+1) for m==0,
    % and Schmidt*sqrt(2l+1) to 4 or 2. Note that for the integration of
    % the harmonics you get either 2pi for m==0 or pi+pi for cosine and
    % sine squared (cross terms drop out). This makes the fully
    % normalized spherical harmonics the only ones that consistently give
    % 1 for the normalization of the spherical harmonics.
    % Note this is using the cosines only; the "spherical harmonics" are
    % actually only semi-normalized.
    % Test normalization as follows (using inaccurate Simpson's rule):


            % Find the cosine and sine coefficients for this degree                         
            % For the Blm
            bcolm=shcos(blmcosi,l);
            bsilm=shsin(blmcosi,l);
            % For the Clm
            ccolm=shcos(clmcosi,l);
            csilm=shsin(clmcosi,l);
                    
            mmult=repmat(m(:),1,length(phi));
        
            % The Blm expansion preparation
            % The derivative by phi part
            % For the negative m, the minus from inside sin and the minus
            % from mmult cancel out
            dphfacBlm=repmat(bcolm,1,nlon).*(-sin(mphi)).*mmult+...
                      repmat(bsilm,1,nlon).*  cos(mphi) .*mmult;

            % The derivative by theta part
            % The minus within cos doesn't affect cos
            dthfacBlm=repmat(bcolm,1,nlon).*cos(mphi)+...
                      repmat(bsilm,1,nlon).*sin(mphi);
              
                  
            % The Clm expansion preparation
            % The derivative by phi part  
            % For the negative m, the minus from inside sin and the minus
            % from mmult cancel out
            dphfacClm=repmat(ccolm,1,nlon).*(-sin(mphi)).*mmult+...
                      repmat(csilm,1,nlon).*  cos(mphi) .*mmult;
              
            % The derivative by theta part  
            % The minus within cos doesn't affect cos
            dthfacClm=repmat(ccolm,1,nlon).*cos(mphi)+...
                      repmat(csilm,1,nlon).*sin(mphi);      
   
            % Sum over all orders and (through loop) over all degrees
            if length(degres)==length(c11cmn)
                
                % The expansion of Blm
                % Theta component
                expaBlmth= sum(dthplm.*dthfacBlm,1)';
                
                % Phi component
                expaBlmph= sum(dphplm.*dphfacBlm,1)';
                
                % The expansion of Clm
                % Theta component
                expaClmth= sum(dphplm.*dphfacClm,1)';
                
                % Phi component
                expaClmph=-sum(dthplm.*dthfacClm,1)';
                                                               
                % Or diag(plm'*fac1) if you will
                
            elseif length(degres)==1 && length(c11cmn)==4
                
                % The expansion of Blm
                % Theta component             
                expaBlmth= dthplm'*dthfacBlm;
                
                % Phi component
                expaBlmph= dphplm'*dphfacBlm;
                
                % The expansion of Clm
                % Theta component
                expaClmth= dphplm'*dphfacClm;
                
                % Phi component
                expaClmph=-dthplm'*dthfacClm;
                
            end

            % The phi component of Blm and Clm combined
            r(:,:,1)=r(:,:,1)+expaBlmph+expaClmph;
            
            % The theta component of Blm and Clm combined 
            r(:,:,2)=r(:,:,2)+expaBlmth+expaClmth;
            
        end
    end
        
    lon=phi*180/pi;
    lat=90-theta*180/pi;
  
    % Prepare output
    vars={r,lon,lat};
    varargout=vars(1:nargout);    
    
end
    
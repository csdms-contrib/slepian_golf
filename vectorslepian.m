function varargout=vectorslepian(Lmax,dom,comp,index,res,c11cmn,C,V,rotb)
% [data,lat,lon,C,V,blmcosi,clmcosi]=vectorslepian(Lmax,dom,comp,index,res,
% c11cmn,C,V,rotb)
%
% Calculates the index best bandlimited vector Slepian function for the
% given domain and maximum bandwidth.
%
% INPUT:
% 
% Lmax      Bandwidth
% dom       The region 'england', 'eurasia', 'australia','greenland', 
%           'africa', 'samerica', 'amazon', 'orinoco', 'gpsnamerica',
%           'antarctica', 'alloceans', 'namerica' [default]
% comp      which component should be calculated 'radial' or 'tangential' 
%           [default: 'tangential']
% index     Which Slepian function should be calculated (the index best). 
%           Default=1 (the best).  
% res       The resolution of the calculated function, or if region + 
%           tangential two different resolutions for the quiver plot
% c11cmn    Corner nodes of lon/lat grid [default: 0 90 360 -90] or in the 
%           polar cap case, the longitude range
% C         The precalculated sorted eigenvector matrix   
% V         The precalculated sorted eigenvalue matrix 
% rotb      Should the theta values for the plot be shifted by 90 degrees
%           (for example to plot antarctica actually at the right position)
%
% OUTPUT: 
%
% data       For radial functions the two dimensional array 
%                with values for longitude and latitude 
%            for tangential functions the three dimensional 
%                array, where data(:,:,1) contains the values of 
%                the phi component and data(:,:,2) of the theta  
%                component. 
%                If res is two dimensional: Different
%                resolutions for the abs value and the quiver plot
% lat       The latitude values (cell for two dim res)
% lon       The latitude values (cell for two dim res)
% C         The sorted eigenvector matrix (to use next time)
% V         The sorted eigenvalue matrix (to use next time)
% blmcosi   Either blmcosi (tangential) or lmcosi (radial), the
%           spherical harmonics coefficients for the sin and cos including 
%           their l and m values for Blm (tangential) or Plm (radial)
% clmcosi   coefficients for Clm (tangential) or [] (radial)
%
% EXAMPLE:
%
% vectorslepian('demo1') Plot the 12 best tangential vector Slepians for
%                        Australia
%
% See also BLMCLM2XYZ, KERNELB, KERNELBP, CAPVECTORSLEPIAN
%
% Last modified by plattner-at-alumni.ethz.ch, 02/05/2013

defval('Lmax',18)
defval('dom','namerica')
defval('comp','tangential')
defval('index',1)
defval('res',1)
defval('C',[])
defval('V',[])
defval('ngl',1000)
defval('c11cmn',[0 90 360 -90]);
defval('rotb',0)

if ~isstr(Lmax)

blmcosi=[];
clmcosi=[];

% The general continent case   

if(isempty(C)||isempty(V))
    % Need to calculate or load the kernel                
    if strcmp(comp,'tangential') 
        % The tangential case
        try
            K=kernelbp(Lmax,dom,[],ngl,rotb);   
        catch
            K=kernelb(Lmax,dom,[],ngl,rotb);  
        end
        filoc=fullfile(getenv('IFILES'),'TANSLEPIANS'); 
    elseif strcmp(comp,'radial') 
        % The radial case
        K=kernelcp(Lmax,dom,[],ngl,rotb);   
        filoc=fullfile(getenv('IFILES'),'RADSLEPIANS'); 
    else
        error('Choose either radial or tangential component')
    end  
    if rotb==1
        fnpl=sprintf('%s/SLEP-%s-%i-1_eig.mat',filoc,dom,Lmax);
    else
        fnpl=sprintf('%s/SLEP-%s-%i_eig.mat',filoc,dom,Lmax);
    end
    if ~(exist(fnpl,'file')==2)        
        [C,V]=eig(K);    
        [V,isrt]=sort(sum(V,1),'descend');
        C=C(:,isrt(1:length(K)));
        save(fnpl,'C','V')
    else
        load(fnpl)
        disp(sprintf('%s loaded by VECTORSLEPIAN',fnpl))
    end
end
% Now we have the sorted eigenvectors and eigenvalues
[demsz,delsz,mz,lmc,mzin]=addmon(Lmax);  

% For the calculation of the data, the two cases must be treated
% separately
if strcmp(comp,'tangential')
    % First, only take those coefficients that belong to the Blm
    CBlm=C(1:size(C,1)/2,:);
    % And those that only belong to the Clm
    CClm=C((size(C,1)/2+1):end,:);    
    % The tangential field does not have an l=0 coefficient. In order to
    % still use the already available functions, proceed as follows
    % First put in dummy values for l=0
    CBlm=[NaN(1,size(C,2));CBlm];    
    CBlm=reshape(insert(CBlm(:,index),0,mzin),2,length(demsz))';   
    CClm=[NaN(1,size(C,2));CClm];    
    CClm=reshape(insert(CClm(:,index),0,mzin),2,length(demsz))';        

    % And now Remove again the l=0 part because it is purely radial
    dems=demsz(2:end);
    dels=delsz(2:end);
    CBlm=CBlm(2:end,:);
    CClm=CClm(2:end,:);    
    blmcosi=[dels dems CBlm];
    clmcosi=[dels dems CClm];
    % Now construct the Slepians from the coefficients
    % on a dense grid 
    if length(res)>1
        [data{1},lon{1},lat{1}]=blmclm2xyz(blmcosi,clmcosi,res(1),c11cmn);
        % and on a less dense grid to show the vectors
        [data{2},lon{2},lat{2}]=blmclm2xyz(blmcosi,clmcosi,res(2),c11cmn);  
    else
        [data,lon,lat]=blmclm2xyz(blmcosi,clmcosi,res,c11cmn);

    end

elseif strcmp(comp,'radial')
    Clm=reshape(insert(C(:,index),0,mzin),2,length(demsz))'; 
    lmcosi=[delsz demsz Clm];
    [data,lon,lat]=plm2xyz(lmcosi,res,c11cmn); 
    blmcosi=lmcosi;

else
    error('Choose either radial or tangential component')         
end
    
    
varns={data,lat,lon,C,V,blmcosi,clmcosi};
varargout=varns(1:nargout);


elseif strcmp(Lmax,'demo1')
    % Plot the 12 best tangential vector slepians for Australia
    % First prepare the variables
    clf;
    fig2print(gcf,'landscape')
    Lmax=20;
    dom='australia'%;
    off=10;
    res=[0.2 3];%[0.1 4];
    C=[]; V=[];
    comp='tangential';
    XY=eval(sprintf('%s(10)',dom));
    XY(:,1)=XY(:,1);
    range=[min(XY(:,1))-off max(XY(:,1))+off min(XY(:,2))-off ...
        max(XY(:,2))+off];
    c11cmn=[range(1) range(4) range(2) range(3)];%upper left, lower right
    fozo=10;               
    disp(sprintf('Shannon is %i',round(2*((Lmax+1)^2-1)*spharea(dom))));
    
    [ah,ha,H]=krijetem(subnum(4,3));
    for index=1:12      
        axes(ah(index))
        hold on 
        % Calculate the index-best Slepain. If the eigenvectors and
        % eigenvalues of the Kernel are already calculated, we can reuse
        % them by entering them as input. If they are not yet calculated or
        % loaded, the function vectorslepian returns them. 
        [data,lat,lon,C,V]=vectorslepian(Lmax,dom,comp,...
            index,res,c11cmn,C,V);
        % First the absolute value plot
        [absdata,dmax]=preparetanplot(data{1});
        %disp(sprintf('Index %d, max abs value = %g',index,...
        %    max(max(abs(data{1})))))
        %ka=absdata/max(max(abs(absdata)));
        % Plot
        imagefnan([range(1) range(4)],[range(2) range(3)],-absdata,...
            'kelicol',[-dmax dmax],[],1,100) 
        % Then the arrows
        %quiver(lon{2}(2:end-1),lat{2}(2:end-1),...
        %    data{2}(2:end-1,2:end-1,1),data{2}(2:end-1,2:end-1,2),'k');
        quiverimage(data{2},lon{2},lat{2},0.05,1);
        % Then the continent
        plot(XY(:,1),XY(:,2),'k')
        axis equal
        axis image
        set(ah(index),'FontSize',fozo-2)
        % Writing out the eigenvalues
        try
            boxtex('ll',ah(index),sprintf('%s = %1.6f','\lambda',...
                V(index)),fozo);   
        end
        hold off
        set(ah(index),'xgrid','off','ygrid','off')
        set(gca,'box','on')
    end
    
    o1=serre(ah(1:3),1.5,'across');
    o2=serre(ah(4:6),1.5,'across');
    o3=serre(ah(7:9),1.5,'across');
    o4=serre(ah(10:12),1.5,'across');
    
    
    
    serre(ha(1:4),[],'down')
    serre(ha(5:8),[],'down')
    serre(ha(9:12),[],'down')
    longticks(ah)
  
    set(ah,'xtick',[100:20:160],'xticklabel',[100:20:160],...
         'ytick',[-50:15:0],'yticklabel',[-50:15:0])
    nolabels(ah(1:9),1); nolabels(ha(5:12),2); deggies(ah)
     
    
    fig2print(gcf,'landscape')
    figna=figdisp(comp,sprintf('%s_%i',dom,Lmax));
    
end


end

function varargout=preparetanplot(data)
% [dat,dmax]=PREPARETANPLOT(data)
%
% Returns the intensity (if data is a struct, then the intensity of 
% data{1})
%
% INPUT:
%
% data  data(:,:,1) is a 2d array containing the values of the vector 
%       field in longitudinal direction,(phi) 
%       data(:,:,2) is a 2d array containing the values of the vector 
%       field in latitudinal direction (theta)
%       OR:
%       data is a struct, where data{1} as mentioned before
% 
% OUTPUT:
%
% dat   2d array containing the intensity values of the vectors (the 
%       negative to make the plotting colors in kelicol red). First 
%       dimension is lat, second is lon)
% dmax  Maximum intensity value
%
% See also KELICOL
%
% Last modified by plattner-at-alumni.ethz.ch, 02/27/2012

if iscell(data)
     absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
else
     absdata=sqrt(data(:,:,1).^2+data(:,:,2).^2);
end

dmax=max(max(absdata)); thresh=100;  
dat=-absdata;
dat(abs(dat)<dmax/thresh)=0;

varns={dat,dmax};
varargout=varns(1:nargout);

end


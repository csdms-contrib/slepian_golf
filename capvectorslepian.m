function varargout=capvectorslepian(Lmax,TH,m,index,res,C,V,Vtot,c11cmn,outdat)
% [data,lat,lon,C,V,Vtot,blmcosi,clmcosi]=capvectorslepian(Lmax,TH,m,index,
% res,C,V,Vtot,c11cmn,outdat)
%
% Calculates the index best bandlimited tangential vector slepian function
% for the polar cap with opening angle TH in degrees. 
%
% If a value for m is provided, the index best tangential Slepian for all 
% 0<l<=m and the given m is calculated.
%
% If no value for m is given, the generally index best tangential Slepian  
% is calculated
%
% INPUT:
% 
% Lmax      Bandwidth
% TH        Opening angle of the polar cap
% m         Specific order, or empty if index best Slepian for all orders
% index     if m is empty, return the index-th best tangential vector 
%           Slepian over all orders. If m is given, return the index-best 
%           tangential vector Slepian of order m. If index is negative and 
%           m is given return the abs(index)-best Slepian for order -abs(m)       
% res       The resolution(s) for the data plot. If given as a vector 
%           [res1 res2] then data is going to be a cell of two data
%           (usefull for example in a quiver plot)
% C         The precalculated matrix (or cell, if m is empty) of 
%           eigenvectors for all indices for m (or for all m if m is empty)      
% V         The precalculated vector (or cell, if m is empty) of 
%           eigenvalues for all indices for m (or for all m, if m is empty)
% Vtot      The precalculated sorted lis t of all eigenvalues together with
%           their m and the rank within one single m.
% c11cmn    Corner nodes of lon/lat grid [default: 0 90-sqrt(eps) 360 -90]
%           OR "lon": a column vector with longitudes [degrees]
% outdat    Calculate the data (1)? or just C, V and Vtot (0) [default: 1] 
%
% OUTPUT:
%
% data      the three dimensional data array, data(:,:,1) contains the phi
%           component and data(:,:,2) the theta component. If res is two-
%           dimensional, then data{1} is the data for the first res entry
%           and data{2} for the second res entry.
% lat       The latitude points for the data array
% lon       The longitude points for the data array
% C         The matrix (given m) or cell of matrices (m empty) of the
%           eigenvectors (vector spherical harmonics coefficients)
% V         The vector (given m) or cell of vectors (m empty) of the
%           eigenvalues (spatial focussing ratio)
% Vtot      If m empty, the complete sorted list of all eigenvalues for all
%           m, including their m and rank within a single m. Format is
%           [eigenvalue  m  rank_within_m]
%
% EXAMPLE:
%
% capvectorslepian('demo1') Plots a random polar cap tangential Slepian
%                           field
%
% capvectorslepian('demo2') Calculate and plot the 32 best concentrated 
%                           absolute values for the tangential vector 
%                           Slepians for the polar cap
%
% capvectorslepian('demo3') Plots a tangential Vector Slepian for m and -m
%
% See also VECTORSLEPIAN, KERNELBM, BLMCLM2XYZ, KERNELTANCAPM, KERNELB,
% KERNELBP
%
% Last modified by plattner-at-alumni.ethz.ch, 08/02/2012
% Small change on 06/21/2016

defval('Lmax',18)
defval('TH',30)
defval('m',[])
defval('index',1)
defval('res',1)
defval('C',[])
defval('V',[])
defval('Vtot',[])
defval('c11cmn',[0 90 360 -90]);
defval('outdat',1)

dirname=fullfile(getenv('IFILES'),'TANCAP');

if ~isstr(Lmax)

if isempty(m)    
     % First figure out if already calculated
     fnpl=fullfile(dirname,sprintf('TANCAP-%i-%i.mat',TH,Lmax));
     if exist(fnpl,'file')==2 
        load(fnpl)
        disp(sprintf('%s loaded by CAPVECTORSLEPIAN',fnpl))
     else
     % Calculate
    if(isempty(C)||isempty(V)) 
        % Solve Slepian's problem for all m=0 to l. Put each list of
        % eigenvalues and matrix of eigenvectors into a cell entry. Make a
        % complete list of the eigenvalues with their according m and
        % internal index (Vtot) and sort it.
        m=0;
        Mm=kerneltancapm(TH,Lmax,m);
        [Cm,Vm]=eig(Mm);    
        [Vm,isrt]=sort(sum(Vm,1),'descend');
        Cm=Cm(:,isrt(1:length(Mm)));
        C{1}=Cm;
        V{1}=Vm;     
        Vtot=[Vm' m*ones(size(Vm,2),1) (1:size(Vm,2))'];
        
        
        for m=1:Lmax
            [Mm,Mmm]=kerneltancapm(TH,Lmax,m);
            [Cm,Vm]=eig(Mm);
            [Cmm,Vmm]=eig(Mmm); 
            [Vm,isrt]=sort(sum(Vm,1),'descend');          
            Cm=Cm(:,isrt(1:length(Mm)));
            [Vmm,isrt]=sort(sum(Vmm,1),'descend');
            Cmm=Cmm(:,isrt(1:length(Mmm)));
            C{2*m}=Cm;
            C{2*m+1}=Cmm; 
            V{2*m}=Vm;  
            V{2*m+1}=Vmm; 
            Vtot=[Vtot;Vm'   m*ones(size(Vm ,2),1) (1:size(Vm ,2))'];
            Vtot=[Vtot;Vmm' -m*ones(size(Vmm,2),1) (1:size(Vmm,2))'];
        end
        % and now sort Vtot
        Vtot=sortrows(Vtot,1);
        Vtot=flipud(Vtot);
    end
    save(fnpl,'C','V','Vtot');
     end
    if outdat
    % The next part is the same as below, but for struct V and struct C
    % Now we have the sorted eigenvectors and eigenvalues
    % Set up a vector such that the entries are distributed to the right m
    % only. Do it for a vector that includes l=0 and then remove the first 
    % element.
    [demsz,delsz,mz,lmc,mzin]=addmon(Lmax); 
    blmcosi=[delsz demsz zeros(length(delsz),1) zeros(length(delsz),1)];
    clmcosi=[delsz demsz zeros(length(delsz),1) zeros(length(delsz),1)];
    % Now figure out which m the index best Slepian has.
    m=Vtot(index,2);
    am=abs(m);
    sg=ceil((1-sign(m))/2); % floor because of m=0
    Cm=C{2*am+sg};
    Vm=V{2*am+sg};
    mindex=Vtot(index,3);  
    % mz contains the positions of the m=0. Therefore mz+ma contains the
    % positions of the m
    mpos=mz+am; 
    % Where the entries are for the specific m. Numbers do not make sense
    % for l<m. Therefore we cut them out in the next line.
    pos=mpos(am+1:end);
    % In the m=0 case we need to take out the l=0 case
    if(m==0)
        pos=pos(2:end);
    end
    if(m>0)
        blmcosi(pos,4)=Cm(1:length(Cm)/2,mindex);
        clmcosi(pos,3)=Cm(length(Cm)/2+1:end,mindex);
    elseif(m<0)
        blmcosi(pos,3)=Cm(1:length(Cm)/2,mindex);
        clmcosi(pos,4)=Cm(length(Cm)/2+1:end,mindex);
    else
        blmcosi(pos,3)=Cm(1:length(Cm)/2,mindex);
        clmcosi(pos,3)=Cm(length(Cm)/2+1:end,mindex);
    end
    % And take away the l=0 part
    blmcosi=blmcosi(2:end,:);
    clmcosi=clmcosi(2:end,:); 
    % Now construct the Slepians from the coefficients
    % on a dense grid 
    if length(res)>1
        [data{1},lon{1},lat{1}]=blmclm2xyz(blmcosi,clmcosi,res(1),c11cmn);
        % and on a less dense grid to show the vectors
        [data{2},lon{2},lat{2}]=blmclm2xyz(blmcosi,clmcosi,res(2),c11cmn);  
    else
        [data,lon,lat]=blmclm2xyz(blmcosi,clmcosi,res,c11cmn);
    end
    else 
        data=[]; blmcosi=[]; clmcosi=[]; lat=[]; lon=[];
    end
else
    % First figure out if already calculated
     fnpl=fullfile(dirname,sprintf('TANCAP-%i-%i-%i.mat',TH,Lmax,m));
     if exist(fnpl,'file')==2 
        load(fnpl)
        disp(sprintf('%s loaded by CAPVECTORSLEPIAN',fnpl))
     else
     % Calculate
    am=abs(m);
    if(m<0&&index>0)
        sig=sign(m);
    else
        sig=sign(index);
    end
    index=abs(index);
    m=sig*am;
    if(index > Lmax-am+1)
        error(sprintf('index must be smaller or equal to Lmax-m = %d',...
            Lmax-m))
    end
    if(isempty(C)||isempty(V))       
        Mm=kerneltancapm(TH,Lmax,m);
        [C,V]=eig(Mm);
        [V,isrt]=sort(sum(V,1),'descend');
        C=C(:,isrt(1:length(Mm)));
    end
     end
    if outdat
    % Now we have the sorted eigenvectors and eigenvalues
    % Set up a vector such that the entries are distributed to the right m
    % only. Do it for a vector that includes l=0 and then remove the first 
    % element.
    [demsz,delsz,mz,lmc,mzin]=addmon(Lmax); 
    blmcosi=[delsz demsz zeros(length(delsz),1) zeros(length(delsz),1)];
    clmcosi=[delsz demsz zeros(length(delsz),1) zeros(length(delsz),1)];
    % mz contains the positions of the m=0. Therefore mz+ma contains the
    % positions of the m
    am=abs(m);
    mpos=mz+am; 
    % Where the entries are for the specific m. Numbers do not make sense
    % for l<m. Therefore we cut them out in the next line.
    pos=mpos(am+1:end);
    % In the m=0 case we need to take out the l=0 case
    if(m==0)
        pos=pos(2:end);
    end
    if(m>0)
        % plattner 5/9/2017:
        % This is where the positive-m-blm and negative-m-clm need to get 
        % switched. 
        % I think it is because in Plattner et al 2014 (ACHA), the
        % D_(lm,l'm') are only non-zero when m'=-m. at the same time,
        % B_lm,l'm=B_l-m,l'-m and C_lm,l'm=C_l-m,l'-m. So ultimately we are
        % solving for the -m-component of Blm when we are solving for the
        % m-component of Clm and vice-versa.
        blmcosi(pos,4)=C(1:length(C)/2,index);
        clmcosi(pos,3)=C(length(C)/2+1:end,index);
    elseif(m<0)
        blmcosi(pos,3)=C(1:length(C)/2,index);
        clmcosi(pos,4)=C(length(C)/2+1:end,index);
    else 
        blmcosi(pos,3)=C(1:length(C)/2,index);
        clmcosi(pos,3)=C(length(C)/2+1:end,index);
    end 
       
    % And take away the l=0 part
    blmcosi=blmcosi(2:end,:);
    clmcosi=clmcosi(2:end,:); 
    % Now construct the Slepians from the coefficients
    % on a dense grid 
    if length(res)>1
        [data{1},lon{1},lat{1}]=blmclm2xyz(blmcosi,clmcosi,res(1),c11cmn);
        % and on a less dense grid to show the vectors
        [data{2},lon{2},lat{2}]=blmclm2xyz(blmcosi,clmcosi,res(2),c11cmn);  
    else
        [data,lon,lat]=blmclm2xyz(blmcosi,clmcosi,res,c11cmn);
    end
    else 
        data=[]; blmcosi=[]; clmcosi=[]; lat=[]; lon=[];
    end
    Vtot=[];
    save(fnpl,'C','V','Vtot');
     
end

varns={data,lat,lon,C,V,Vtot,blmcosi,clmcosi};
varargout=varns(1:nargout);

elseif strcmp(Lmax,'demo1')
    TH=20+ceil(rand*40);
    Lmax=10+ceil(rand*(20))
    m=floor(rand*(2*Lmax+1))-Lmax
    index=ceil(rand*(Lmax-abs(m)))
    res=[0.5 2];% Choose [0.3 2] as res for a nicer image, [1 2] for faster
    [data,lat,lon,C,V]=capvectorslepian(Lmax,TH,m,index,res);
    [dat,dmax]=preparetanplot(data);
    clf;
    hold on
    [~,ch,~]=plotplm(dat,lon,lat,5,res(1),TH,1,[-dmax dmax]);
    delete(ch);
    quivpolcaps(data{2})
    axis tight
    try
      boxtex('ll',gca,sprintf('%s =%1.3f','\lambda',V(index)),15);
    end
    hold off   
    
elseif strcmp(Lmax,'demo2')       
    Lmax=18;
    TH=40;
    m=[];
    res=1; % To make it nicer, set res=0.2;
    C=[];
    Vtot=[];
    clf;
    fig2print(gcf,'tall')
    [ah,ha,H]=krijetem(subnum(8,4));        
    C=[];V=[];Vtot=[];
    ind=1;
    index=1;
    while ind<=32
        axes(ah(ind));
        [C,V,Vtot]=plotfriedegg(Lmax,TH,index,res,C,V,Vtot);        
        tl(ind)=title(sprintf('%s_{%i} =%1.3f; m = %i',...
           '\lambda',Vtot(index,3),Vtot(index,1),Vtot(index,2)));
        if(Vtot(index,2)>=0)
            ind=ind+1;
        end
       index=index+1;
    end
    seemax(ah,3)
    set(ah,'CameraV',6.25)
    % Last-minute cosmetics
    movev(tl,-0.25)
    figna=figdisp('tangential',sprintf('polarcap_abs_L%i',Lmax));     
    
elseif strcmp(Lmax,'demo3')
    Lmax=18;
    direction=1; % 0=vertical, 1=horizontal
    TH=40;
    m=1;
    res=[0.2 2];%[1 5];%[0.3 2];
    V=[]; C=[]; Vtot=[];
    %clf;
    thresh=100;
    if direction
        [ah,ha,H]=krijetem(subnum(1,2));
    else
        [ah,ha,H]=krijetem(subnum(2,1));
    end
    index=1;
    
    % Blm,Cl-m
    [data,lat,lon,C,V,Vtot]=capvectorslepian(Lmax,TH,m,index,res);
    axes(ah(1))
    [dat,dmax]=preparetanplot(data);%
    dat(abs(dat)<dmax/thresh)=0;
    [~,ch,~]=plotplm(dat,lon,lat,5,res(1),TH,1,[-dmax dmax]);
    delete(ch);
    quivpolcaps(data{2})     
    axis tight 
    tl(1)=title(sprintf('%s_{%i} =%1.3f; m = %d','\lambda',...
        index,V(index),m));                


    % Bl-m,Clm
    m=-m;
    [data,lat,lon,C,V,Vtot]=capvectorslepian(Lmax,TH,m,index,res);
    axes(ah(2))
    [dat,dmax]=preparetanplot(data);
    dat(abs(dat)<dmax/thresh)=0;
    [~,ch,~]=plotplm(dat,lon,lat,5,res(1),TH,1,[-dmax dmax]);
    delete(ch);
    quivpolcaps(data{2})      
    tl(2)=title(sprintf('%s_{%i} =%1.3f; m = %d','\lambda',...
        index,V(index),m));       
    axis tight    
    
    set(tl,'FontS',12)
    
    if direction
        serre(ah(1:2),0.5,'across')
    else
         serre(ah(1:2),0.5,'down')
    end
    if direction
        fig2print(gcf,'landscape')
        figna=figdisp('tangential',sprintf(...
            'polarcap-%i_L-%i_m-%i_horiz',TH,Lmax,abs(m)),[],1,'pdf'); 
    else
        fig2print(gcf,'flandscape')
        figna=figdisp('tangential',sprintf(...
            'polarcap-%i_L-%i_m-%i_vert',TH,Lmax,abs(m)),[],1,'pdf'); 
    end
    
end

%%%%%%%%%%%%%%% For the plotting in demo2 %%%%%%%%%%%%%%%%
function [C,V,Vtot]=plotfriedegg(Lmax,TH,index,res,C,V,Vtot)

[vdata,lat,lon,C,V,Vtot]=capvectorslepian(Lmax,TH,[],index,...
    res,C,V,Vtot);                        
data=sqrt(vdata(:,:,1).^2+vdata(:,:,2).^2);            
cbb=gray(10);
cbb=kelicol;
ispl=-data/max(abs(data(:)));
ispl(abs(ispl)<1/100)=NaN;              
[~,ch,~]=plotplm(ispl,lon,lat,5,res(1),TH,1,[-1 1]);
delete(ch);
colormap(cbb); hold on
ax2(1)=circ(1); ax2(2)=circ(sin(TH*pi/180));
set(ax2(1:2),'LineW',1)
set(ax2(2),'LineS','--')
axis([-1.0100    1.0100   -1.0100    1.0100])
axis off
     
      

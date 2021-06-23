function varargout=vectorspectral(L,dom,Llim,comp,index,res,c11cmn)
% [hdata,coeff,coeffmatrix,V]=vectorspectral(L,dom,Llim,comp,index,res,
% c11cmn)
%
% Given a domain and a bandwidth, this function calculates the index best
% spatially limited, and spectrally optimized functions by truncating the 
% index best spectrally limited and spatially optimized function and then
% calculating it's spectrum (See for example SDW 2006).
%
% INPUT:
%
% L             Bandwidth to which the function shall be optimized
% Llim          Maximum bandwidth up to which the coefficients should be
%               calculated.
% comp          The component: 'radial' or 'tangential'
% index         The Slepian index (the "index"-best function)
% res           The resolution for the spatial truncation
%
% OUTPUT:
% 
% hdata         The 'index'-best spatially truncated spectrally optimized 
%               Slepian (point data over the globe)
% coeff         The coefficients of the Slepian up to Llim
% coeffmatrix   The coefficients arranged in a matrix that is easy to plot
% V             The eigenvalues (spectral optimization ratios)
%
% EXAMPLE: 
%
% vectorspectral('demo1') Spectra of the 12 best spectrally optimized 
%                         Slepians for North America, the radial component
%                                  
% vectorspectral('demo2') Spectra of the 9 best spectrally optimized
%                         Slepians for North America, the tangential
%                         components blm and clm
%
% vectorspectral('demo3') Comparison of the spatially optimized Slepian to
%                         the spectrally optimized Slepian, the radial 
%                         component
%
% vectorspectral('demo4') Same as demo3 the tangential components blm and 
%                         clm
%                                  
% vectorspectral('demo5') Whole world comparison of the spatially optimized 
%                         Slepian to the spectrally optimized Slepian, the
%                         tangential components blm and clm
%
% vectorspectral('demo6') Same as demo5 for the radial component
%
% vectorspectral('demo7') Comparing the predicted focussing ratio 
%                         (eigenvalue) to the numerically calculated
%                         spatial and spectral energy focussing and to
%                         the ratio of the spherical harmonics coefficients
%                         the radial component
% 
% vectorspectral('demo8') Same as demo7 for the tangential components blm
%                         clm
%
% See also VECTORSLEPIAN, ERRORINREGION
%
% Last modified by plattner-at-alumni.ethz.ch, 02/28/2012


defval('dom','namerica')
defval('L',10)
defval('Llim',20)
defval('comp','radial')
defval('index',1)
defval('res',1)
defval('c11cmn',[0 90-sqrt(eps) 360 -90])

if ~isstr(L)
% First we calculate the spectrally limited and spatially optimized Slepian
% function g_index
[data,lat,lon,C,V]=vectorslepian(L,dom,comp,index,res,c11cmn);
% Then the index best spatially limited and spectrally optimized Slepian is
% equal to the spatially limited g_index
if length(res)==1
    hdata=errorinregion(-data,[],dom,c11cmn);
elseif length(res)==2
    hdata{1}=errorinregion(-data{1},[],dom,c11cmn);
    hdata{2}=errorinregion(-data{2},[],dom,c11cmn);
else
    error('Too many entries in res')
end

% If we want to know the spectrum of h_index, we need to decompose hdata
% into spherical harmonics
[demsz,delsz,mz,lmc,mzin,mzo]=addmon(Llim);
if strcmp(comp,'radial')
    lmcosi=xyz2plm(hdata,Llim);
    [coeffmatrix,coeff]=intomatrix(lmcosi,Llim);    
elseif strcmp(comp,'tangential')
    if length(res)==1
        [blmcosi, clmcosi]=xyz2blmclm(hdata,Llim);
    elseif length(res)==2
        % Assuming that res(1)>res(2)
        [blmcosi, clmcosi]=xyz2blmclm(hdata{1},Llim);
    else
        error('Too many entries in res')
    end  
    [bcoeffmatrix,bcoeff]=intomatrix(blmcosi,Llim,0);
    [ccoeffmatrix,ccoeff]=intomatrix(clmcosi,Llim,0);    
    coeff{1}=bcoeff;
    coeff{2}=ccoeff;
    coeffmatrix{1}=bcoeffmatrix;
    coeffmatrix{2}=ccoeffmatrix;    
else
    error('Choose either "tangential" or "radial" as component')
end

varns={hdata,coeff,coeffmatrix,V};
varargout=varns(1:nargout);

elseif strcmp(L,'demo1')
    L=20;
    Llim=50;
    dom='namerica';
    comp='radial';
    res=1;    
    [ah,ha,H]=krijetem(subnum(4,3));    
    for index=1:12
        [hdata,coeff,coeffmatrix,V]=vectorspectral(L,dom,Llim,...
        comp,index,res);
        axes(ah(index))
        imagefnan([-Llim,0],[Llim,Llim],coeffmatrix,[],[],[],1,100);
        hold on
        plot([-Llim Llim],[L+0.5 L+0.5],'k')
	try
          boxtex('ll',ah(index),sprintf('%s =%1.3f','\lambda',V(index)),5);
	end
        hold off
    end
    serre(ah(1:3),2/3,'across')
    serre(ah(4:6),2/3,'across')
    serre(ah(7:9),2/3,'across')
    serre(ah(10:12),2/3,'across')
    serre(ha(1:4),2/3,'down')
    serre(ha(5:8),2/3,'down')   
    serre(ha(9:12),2/3,'down')
    nolabels(ah(1:9),1); nolabels(ha(5:12),2);
    figdisp('spectral',sprintf('%s_%s_%i',comp,dom,L))
    
elseif strcmp(L,'demo2')    
    L=20;
    Llim=50;
    dom='australia';
    comp='tangential';
    res=1;  
    range=[0 360 -90 90-sqrt(eps)];%[190 550 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];     
    [ah,ha,H]=krijetem(subnum(3,3));    
    figure
    [bh,hb,B]=krijetem(subnum(3,3));
    for index=1:9
        [hdata,coeff,coeffmatrix,V]=vectorspectral(L,dom,Llim,...
        comp,index,res,c11cmn);
        % Plot blm coefficients 
        axes(ah(index))
        imagefnan([-Llim,0],[Llim,Llim],coeffmatrix{1},[],[],[],1,100);
        hold on
        xlabel('m')
        ylabel('L')
        plot([-Llim Llim],[L+0.5 L+0.5],'k')
  try
          boxtex('ll',ah(index),sprintf('%s =%1.3f','\lambda',V(index)),5);
	end
	hold off
        % Plot clm coefficients 
        axes(bh(index))
        imagefnan([-Llim,0],[Llim,Llim],coeffmatrix{2},[],[],[],1,100);
        hold on
        xlabel('m')
        ylabel('L')
        plot([-Llim Llim],[L+0.5 L+0.5],'k')
	try
          boxtex('ll',bh(index),sprintf('%s =%1.3f','\lambda',V(index)),5);
	end
        hold off
    end
    
elseif strcmp(L,'demo3')  
    index=1;
    L=20;
    if ~isstr(dom)
        L=dom
        index=Llim
    end  
    Llim=2*L;%50;
    dom='namerica';
    comp='radial';
    res=1;     
    clf
    XY=eval(sprintf('%s(10)',dom));
    off=10;
    range=[0 360 -90 90];
    c11cmn=[range(1) range(4) range(2) range(3)];%upper left, lower right    
    [ah,ha,H]=krijetem(subnum(2,2));
    % First the spectrally limited spatially optimized Slepian ...
    data=vectorslepian(L,dom,comp,index,res);
    axes(ah(1))
    imagefnan([range(1) range(4)],[range(2) range(3)],data);
    hold on
    plot(XY(:,1),XY(:,2),'k')
    hold off
    % ... with it's spectral content
    lmcosi=xyz2plm(data,Llim);
    [coeffmatrix,coeff]=intomatrix(lmcosi,Llim); 
    axes(ah(2))
    imagefnan([-Llim,0],[Llim,Llim],coeffmatrix);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')
    hold off
    % Then the spatially limited spectrally optimized Slepian ...
    [hdata,coeff,coeffmatrix,V]=vectorspectral(L,dom,Llim,...
        comp,index,res);
    axes(ah(3))
    imagefnan([range(1) range(4)],[range(2) range(3)],hdata);
    hold on
    plot(XY(:,1),XY(:,2),'k')
    hold off
    % ... with it's spectral content
    axes(ah(4))
    imagefnan([-Llim,0],[Llim,Llim],coeffmatrix);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')      
    hold off
    try
      boxtex('ll',ah(1),sprintf('%s =%1.3f','\lambda',V(index)),12); 
      boxtex('ll',ah(4),sprintf('%s =%1.3f','\lambda',V(index)),12);
    end
    
    serre(ha(1:2),3/2,'down')
    serre(ha(3:4),3/2,'down')      
    nolabels(ah(1:2),1); 
    figdisp('spectral',sprintf('%s_%s_%i_%i',comp,dom,L,index))
    
    
    
    
elseif strcmp(L,'demo4')  
    index=20;
    L=20;
    if ~isstr(dom)
        L=dom
        index=Llim
    end  
    Llim=2*L;
    dom='samerica';
    comp='tangential';
    res=[.1 5];       
    clf
    XY=eval(sprintf('%s(10)',dom));
    off=10;
    range=[0 360 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];%upper left, lower right    
    [ah,ha,H]=krijetem(subnum(2,3));
    % First the spectrally limited spatially optimized Slepian ...
    [data,lat,lon]=vectorslepian(L,dom,comp,index,res,c11cmn); 
    absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
    axes(ah(1))
    imagefnan([range(1) range(4)],[range(2) range(3)],-absdata);
    hold on
    quiverimage(data{2},lon{2},lat{2},0.05)
    plot(XY(:,1),XY(:,2),'k')
    hold off
    % ... with it's spectral content
    [blmcosi clmcosi]=xyz2blmclm(data{1},Llim);
    % Get rid of machine precision error
    blmcosi(abs(blmcosi(:,3))<sqrt(eps),3)=0;
    blmcosi(abs(blmcosi(:,4))<sqrt(eps),4)=0;
    clmcosi(abs(clmcosi(:,3))<sqrt(eps),3)=0;
    clmcosi(abs(clmcosi(:,4))<sqrt(eps),4)=0;
    bcoeffmatrix=intomatrix(blmcosi,Llim,0); 
    ccoeffmatrix=intomatrix(clmcosi,Llim,0); 
    % blm
    axes(ah(2))
    imagefnan([-Llim,0],[Llim,Llim],bcoeffmatrix);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')  
    hold off
    % clm
    axes(ah(3))
    imagefnan([-Llim,0],[Llim,Llim],ccoeffmatrix);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')  
    hold off
    % Then the spatially limited spectrally optimized Slepian ...
    [hdata,coeff,coeffmatrix,V]=vectorspectral(L,dom,Llim,...
        comp,index,res);
    abshdata=sqrt(hdata{1}(:,:,1).^2+hdata{1}(:,:,2).^2);
    axes(ah(4))
    imagefnan([range(1) range(4)],[range(2) range(3)],-abshdata);
    hold on
    plot(XY(:,1),XY(:,2),'k')
    quiverimage(hdata{2},lon{2},lat{2},0.05)
    hold off
    % ... with it's spectral content
    % blm
    axes(ah(5))
    imagefnan([-Llim,0],[Llim,Llim],coeffmatrix{1});
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')      
    hold off
    % clm
    axes(ah(6))
    imagefnan([-Llim,0],[Llim,Llim],coeffmatrix{2});
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')      
    hold off
    try
      boxtex('ll',ah(1),sprintf('%s =%1.3f','\lambda',V(index)),7); 
      boxtex('ll',ah(5),sprintf('%s =%1.3f','\lambda',V(index)),7);      
      boxtex('ll',ah(6),sprintf('%s =%1.3f','\lambda',V(index)),7);
    end
    serre(ha(1:2),4.5/2,'down')
    serre(ha(3:4),4.5/2,'down')  
    serre(ha(5:6),4.5/2,'down')
    serre(ah(1:2),3/3,'across')
    serre(ha(2:3),3/3,'across')
    serre(ah(2:3),2/3,'across')
    serre(ah(5:6),2/3,'across')
    nolabels(ah(3),2); 
    nolabels(ha(6),2); 
    nolabels(ah(1:3),1); 
    figdisp('spectral',sprintf('%s_%s_%i_%i',comp,dom,L,index))    
    
elseif strcmp(L,'demo5')  
    % The world tangential
    doms={'africa', 'eurasia', 'namerica', 'australia', 'greenland', ...
      'samerica'};
    dom='world';
    index=160;
    L=20; 
    Llim=2*L;
    comp='tangential';
    res=[0.1 10];%[.5 10];   
    fs=14;
    del=0.1;
    quiverstretch=0;
    fact=3; 
    
    range=[190 550 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)]; %upper left, lower right    
            
    % First get the Kernel matrices for all eigenvalues
    % Check if already calculated for the world
    filoc=fullfile(getenv('IFILES'),'TANSLEPIANS'); 
    fnpl=sprintf('%s/SLEP-%s-%i_eig.mat',filoc,dom,L);
    if ~(exist(fnpl,'file')==2)   
        K=kernelbp(L,doms{1});
        for d=2:length(doms)                   
            K=K+kernelbp(L,doms{d});             
        end
        [C,V]=eig(K);
        [V,isrt]=sort(sum(V,1),'descend');
        C=C(:,isrt(1:length(K)));
        save(fnpl,'C','V')
    else
       load(fnpl)
       disp(sprintf('%s loaded by VECTORSPECTRAL',fnpl))
    end
    
    clf 
    [ah,ha,H]=krijetem(subnum(3,2));
    % First the spectrally limited spatially optimized Slepian ...
    [data,lat,lon]=vectorslepian(L,dom,comp,index,res,c11cmn,C,V); 
    absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
    axes(ha(1))
    cax=[-max(max(absdata)) max(max(absdata))];
    imagefnan([range(1) range(4)],[range(2) range(3)],-absdata,[],...
        cax,[],[],100);
    hold on   
    quiverimage(fact*data{2},lon{2},lat{2},del,quiverstretch) 
    area=0;
    for d=1:length(doms)
        XY=eval(sprintf('%s(10)',doms{d}));
        area=area+spharea(doms{d});
        if strcmp(doms{d},'australia')
            XY(:,1)=XY(:,1)+360;
        end
        plot(XY(:,1),XY(:,2),'k')
    end
    disp(sprintf('Area of the combined continents is %f',area));
    hold off
    % ... with it's spectral content
    [blmcosi clmcosi]=xyz2blmclm(data{1},Llim);
    % Get rid of machine precision error
    blmcosi(abs(blmcosi(:,3))<sqrt(eps),3)=0;
    blmcosi(abs(blmcosi(:,4))<sqrt(eps),4)=0;
    clmcosi(abs(clmcosi(:,3))<sqrt(eps),3)=0;
    clmcosi(abs(clmcosi(:,4))<sqrt(eps),4)=0;
    bcoeffmatrix=intomatrix(blmcosi,Llim,0); 
    ccoeffmatrix=intomatrix(clmcosi,Llim,0); 
    longticks(ha(1))
    % blm
    axes(ha(2))
    imagefnan([-Llim,0],[Llim,Llim],bcoeffmatrix,[],[],[],1,100);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')  
    hold off
    ylabel('degree l','FontSize',fs)
    longticks(ha(2))
    % clm
    axes(ha(3))
    imagefnan([-Llim,0],[Llim,Llim],ccoeffmatrix,[],[],[],1,100);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')  
    hold off
    % Then the spatially limited spectrally optimized Slepian ...
    % Here I replace the call of vectorslepian because my domain consists 
    % of several domains.
    hdata{1}=errorinregion(-data{1},[],doms{1},c11cmn);
    hdata{2}=errorinregion(-data{2},[],doms{1},c11cmn);
    for d=2:length(doms)
        if strcmp(doms{d},'australia')
            hdata{1}=hdata{1}+errorinregion(-data{1},[],doms{d},...
                c11cmn-[360 0 360 0]);
            hdata{2}=hdata{2}+errorinregion(-data{2},[],doms{d},...
                c11cmn-[360 0 360 0]);
        else
            hdata{1}=hdata{1}+errorinregion(-data{1},[],doms{d},c11cmn);
            hdata{2}=hdata{2}+errorinregion(-data{2},[],doms{d},c11cmn);
        end
    end
    [blmcosi, clmcosi]=xyz2blmclm(hdata{1},Llim);
    [bcoeffmatrix,bcoeff]=intomatrix(blmcosi,Llim,0);
    [ccoeffmatrix,ccoeff]=intomatrix(clmcosi,Llim,0);      
    abshdata=sqrt(hdata{1}(:,:,1).^2+hdata{1}(:,:,2).^2);
    xlabel('order m','FontSize',fs)
    ylabel('degree l','FontSize',fs)
    longticks(ha(3))
    axes(ha(4))
    imagefnan([range(1) range(4)],[range(2) range(3)],-abshdata,[],...
        cax,[],[],100);
    hold on
    for d=1:length(doms)
        XY=eval(sprintf('%s(10)',doms{d}));
        if strcmp(doms{d},'australia')
            XY(:,1)=XY(:,1)+360;
        end
        plot(XY(:,1),XY(:,2),'k')
    end
    quiverimage(fact*hdata{2},lon{2},lat{2},del,quiverstretch)
    hold off
    longticks(ha(4))
    % ... with it's spectral content
    % blm
    axes(ha(5))
    imagefnan([-Llim,0],[Llim,Llim],bcoeffmatrix,[],[],[],1,100);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')      
    hold off
    longticks(ha(5))
    % clm
    axes(ha(6))
    imagefnan([-Llim,0],[Llim,Llim],ccoeffmatrix,[],[],[],1,100);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')      
    hold off
    xlabel('order m','FontSize',fs)  
    try
    [bhan, thhan]=boxtex('ll',ha(2),sprintf('%s','{V}_{lm} '),fs-2,[],1.2); 
    delete(bhan)
    [bhan, thhan]=boxtex('ll',ha(3),sprintf('%s','{W}_{lm} '),fs-2,[],1.2); 
    delete(bhan)
    [bhan, thhan]=boxtex('ll',ha(5),sprintf('%s','{V\prime}_{lm} '),...
        fs-2,[],1.2); 
    delete(bhan)
    [bhan, thhan]=boxtex('ll',ha(6),sprintf('%s','{W\prime}_{lm} '),...
        fs-2,[],1.2);  
    delete(bhan)
    end
    longticks(ha(6))
    try
      boxtex('lr',ha(1),sprintf('%s =%1.3f','\lambda',V(index)),fs); 
      boxtex('lr',ha(5),sprintf('%s =%1.3f','\lambda',V(index)),fs);    
      boxtex('lr',ha(6),sprintf('%s =%1.3f','\lambda',V(index)),fs);    
    end
    set([ah(1:2)],'xtick',[270:90:540],'xtickl',[270 0 90 180],...
	 'ytick',[-90:45:90],'ytickl',[-90:45:90])
    set([ah(3:6)],'xtick',[-40:20:40],'xtickl',[-40:20:40],...
	 'ytick',[0:10:40],'ytickl',[0:10:40])
    
    deggies(ah(1:2))
    serre(ha(1:2),3.5,'down')
    serre(ha(2:3),1.8,'down')
    serre(ha(4:5),3.5,'down')
    serre(ha(5:6),1.8,'down')
    serre(ah(1:2),0.7,'across')
    serre(ah(3:4),0.7,'across')
    serre(ah(5:6),0.7,'across')
    nolabels(ha(4:6),2); 
    nolabels(ah(3:4),1); 
    
    fig2print(gcf,'tall')     
    figdisp('spectral',sprintf('%s_%s_%i_%i_switch',comp,dom,L,index),[],1)
    
elseif strcmp(L,'demo6')  
    % The world radial
    doms={'africa', 'eurasia', 'namerica', 'australia', 'greenland', ...
      'samerica'};
    dom='world';
    index=80;
    L=20; 
    Llim=2*L;
    comp='radial';
    res=.1;   
    fs=14;
    
    range=[190 550 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];%upper left, lower right    
            
    % First get the Kernel matrices for all eigenvalues
    % Check if already calculated for the world
    filoc=fullfile(getenv('IFILES'),'RADSLEPIANS'); 
    fnpl=sprintf('%s/SLEP-%s-%i_eig.mat',filoc,dom,L);
    if ~(exist(fnpl,'file')==2)   
        K=kernelcp(L,doms{1});
        for d=2:length(doms)                   
            K=K+kernelcp(L,doms{d});             
        end
        N=sum(diag(K));
        disp(sprintf('Shannon from sum trace of radial kernel is %f',N));
        [C,V]=eig(K);
        [V,isrt]=sort(sum(V,1),'descend');
        C=C(:,isrt(1:length(K)));
        save(fnpl,'C','V')
    else
       load(fnpl)
       N=sum(V);
       disp(sprintf('Shannon from sum trace of radial kernel is %f',N));
       disp(sprintf('%s loaded by VECTORSLEPIAN',fnpl))
    end

    clf 
    [ah,ha,H]=krijetem(subnum(2,2));
    % First the spectrally limited spatially optimized Slepian ...
    [data,lat,lon]=vectorslepian(L,dom,comp,index,res,c11cmn,C,V); 
    axes(ah(1))
    cax=[-max(max(abs(data))) max(max(abs(data)))];   
    imagefnan([range(1) range(4)],[range(2) range(3)],-data,[],cax,...
        [],[],100);
    hold on   
    for d=1:length(doms)
        XY=eval(sprintf('%s(10)',doms{d}));
        if strcmp(doms{d},'australia')
            XY(:,1)=XY(:,1)+360;
        end
        plot(XY(:,1),XY(:,2),'k');
    end
    hold off
    % ... with it's spectral content
    lmcosi=xyz2plm(data,Llim);
    % Get rid of machine precision error
    lmcosi(abs(lmcosi(:,3))<sqrt(eps),3)=0;
    lmcosi(abs(lmcosi(:,4))<sqrt(eps),4)=0;
    coeffmatrix=intomatrix(lmcosi,Llim); 
    longticks(ah(1))
    axes(ah(3))
    imagefnan([-Llim,0],[Llim,Llim],coeffmatrix,[],[],[],[],100);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')  
    hold off 
    xlabel('order m','FontSize',fs)
    ylabel('degree l','FontSize',fs)
    % Then the spatially limited spectrally optimized Slepian ...
    % Here I replace the call of vectorslepian because my domain consists 
    % of several domains.
    hdata=errorinregion(-data,[],doms{1},c11cmn);
    for d=2:length(doms)
        if strcmp(doms{d},'australia')
           hdata=hdata+errorinregion(-data,[],doms{d},...
               c11cmn-[360 0 360 0]);
        else
           hdata=hdata+errorinregion(-data,[],doms{d},c11cmn); 
        end       
    end
    lmcosi=xyz2plm(hdata,Llim);
    [coeffmatrix,coeff]=intomatrix(lmcosi,Llim);      
    longticks(ah(3))
    axes(ah(2))
    imagefnan([range(1) range(4)],[range(2) range(3)],-hdata,[],...
        cax,[],[],100);
    hold on
    for d=1:length(doms)
        XY=eval(sprintf('%s(10)',doms{d}));
        if strcmp(doms{d},'australia')
            XY(:,1)=XY(:,1)+360;
        end
        plot(XY(:,1),XY(:,2),'k')
    end
    hold off
    longticks(ah(2))
    % ... with it's spectral content
    axes(ah(4))
    imagefnan([-Llim,0],[Llim,Llim],coeffmatrix);
    hold on    
    plot([-Llim Llim],[L+0.5 L+0.5],'k')      
    hold off  
    xlabel('order m','FontSize',fs)
    longticks(ah(4))
    try
    [bhan, thhan]=boxtex('ll',ha(2),sprintf('%s','{U}_{lm} '),fs-2,[],1.2); 
    delete(bhan)   
    [bhan, thhan]=boxtex('ll',ha(4),sprintf('%s','{U\prime}_{lm} '),...
        fs-2,[],1.2); 
    delete(bhan)  
    
    boxtex('lr',ah(1),sprintf('%s =%1.3f','\lambda',V(index)),fs); 
    boxtex('lr',ah(4),sprintf('%s =%1.3f','\lambda',V(index)),fs);     
    end
    serre(ha(1:2),3/2,'down')
    serre(ha(3:4),3/2,'down')  
    set([ah(1:2)],'xtick',[270:90:540],'xtickl',[270 0 90 180],...
	 'ytick',[-90:45:90],'ytickl',[-90:45:90])
    set([ah(3:4)],'xtick',[-40:20:40],'xtickl',[-40:20:40],...
	 'ytick',[0:10:40],'ytickl',[0:10:40])
    serre(ah(1:2),0.7,'across')
    serre(ah(3:4),0.7,'across')
    deggies(ah(1))
    deggies(ah(2))
    nolabels(ah(4),2); 
    nolabels(ah(2),2);

    fig2print(gcf,'portrait')
    figdisp('spectral',sprintf('%s_%s_%i_%i_switch',comp,dom,L,index),[],1)  
    
elseif strcmp(L,'demo7')  
    L=20;
    Llim=100;
    dom='africa';
    comp='radial';
    res=.01;    
    index=20;
    range=[190 550 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];    
    gdata=vectorslepian(L,dom,comp,index,res,c11cmn);
    lmcosi=xyz2plm(gdata,Llim);
    [gcoeffmatrix,gcoeff]=intomatrix(lmcosi,Llim);     
    [hdata,hcoeff,hcoeffmatrix,V]=vectorspectral(L,dom,Llim,...
        comp,index,res,c11cmn);
    % Turn the NaNs (good for plotting) into zeros (good for summing)
    nans=find(isnan(hcoeffmatrix));
    hcoeffmatrix(nans)=0;    
    % The energy norm of the coefficients in total is
    totenergy=sum(sum(hcoeffmatrix.^2));
    % The energy norm of the coefficients the target range is 
    Lenergy=sum(sum(hcoeffmatrix(1:(L+1),:).^2)); 
    % Going until L+1 because indexing starts with 1
    
    disp(sprintf(...
       'Energy ratio L: Theor: %g, Actual: %g (measured until Llim=%d)',...
       V(index),Lenergy/totenergy,Llim))    
    totenergy2=sum(sum(gdata.^2));
    Renergy=sum(sum(hdata.^2));    
    disp(sprintf('Energy ratio R: Theor: %g, Actual: %g',...
        V(index),Renergy/totenergy2))      
    ratio=hcoeffmatrix(1:(L+1),(Llim+1-L):(Llim+1+L))...
        ./gcoeffmatrix(1:(L+1),(Llim+1-L):(Llim+1+L));   
    cax=V(index)*[0.75 1.25];
    imagesc([-L,L],[0,L],ratio,cax)
    axis xy
    colorbar
    title(sprintf('Ratio h_{lm}/g_{lm}, should be %g',V(index)))
    xlabel('m')
    ylabel('l')    
    nans=find(isnan(ratio));
    ratio(nans)=0;
    avgratio=sum(sum(abs(ratio)))/(L+1)/(L+1);
    disp(sprintf('Average ratio = %g',avgratio));   
    figdisp('spectral_ratio',sprintf('%s_%s_%i_%i',comp,dom,L,index))
    
    
    elseif strcmp(L,'demo8')  
    L=20;
    Llim=100;
    dom='africa';
    comp='tangential';
    res=.1;    
    index=100;
    range=[190 550 -90 90-sqrt(eps)];
    c11cmn=[range(1) range(4) range(2) range(3)];    
    gdata=vectorslepian(L,dom,comp,index,res,c11cmn);
    [blmcosi,clmcosi]=xyz2blmclm(gdata,Llim);
    [gcoeffmatrix_blm,gcoeff]=intomatrix(blmcosi,Llim,0);  
    [gcoeffmatrix_clm,gcoeff]=intomatrix(clmcosi,Llim,0);  
    [hdata,hcoeff,hcoeffmatrix,V]=vectorspectral(L,dom,Llim,...
        comp,index,res,c11cmn);
    % Turn the NaNs (good for plotting) into zeros (good for summing)
    nans=find(isnan(hcoeffmatrix{1}));
    hcoeffmatrix{1}(nans)=0;    
    nans=find(isnan(hcoeffmatrix{2}));
    hcoeffmatrix{2}(nans)=0;   
    % The energy norm of the coefficients in total is
    totenergy_blm=sum(sum(hcoeffmatrix{1}.^2));
    totenergy_clm=sum(sum(hcoeffmatrix{2}.^2));
    % The energy norm of the coefficients the target range is 
    Lenergy_blm=sum(sum(hcoeffmatrix{1}(1:(L+1),:).^2)); 
    Lenergy_clm=sum(sum(hcoeffmatrix{2}(1:(L+1),:).^2)); 
    % Going until L+1 because indexing starts with 1
    
    disp(sprintf(...
   'Energy ratio L Blm: Theor: %g, Actual: %g (measured until Llim=%d)',...
   V(index),Lenergy_blm/totenergy_blm,Llim))    
    disp(sprintf(...
   'Energy ratio L Clm: Theor: %g, Actual: %g (measured until Llim=%d)',...
   V(index),Lenergy_clm/totenergy_clm,Llim))      
    totenergy2=sum(sum(gdata(:,:,1).^2+gdata(:,:,2).^2));
    Renergy=sum(sum(hdata(:,:,1).^2+hdata(:,:,2).^2));    
    disp(sprintf('Energy ratio R: Theor: %g, Actual: %g',...
        V(index),Renergy/totenergy2))      
    ratio_blm=hcoeffmatrix{1}(1:(L+1),(Llim+1-L):(Llim+1+L))...
        ./gcoeffmatrix_blm(1:(L+1),(Llim+1-L):(Llim+1+L));   
    ratio_clm=hcoeffmatrix{2}(1:(L+1),(Llim+1-L):(Llim+1+L))...
        ./gcoeffmatrix_clm(1:(L+1),(Llim+1-L):(Llim+1+L));   
    %ratio=setnans(ratio,0);
    [ah,ha,H]=krijetem(subnum(1,3));
    
    axes(ah(1))
    cax=V(index)*[0.75 1.25];
    imagesc([-L,L],[0,L],ratio_blm,cax)
    axis xy
    colorbar
    title(sprintf('Ratio B_{lm}: h_{lm}/g_{lm}, should be %g',V(index)))
    xlabel('m')
    ylabel('l')    
    
    axes(ah(2))
    cax=V(index)*[0.75 1.25];
    imagesc([-L,L],[0,L],ratio_clm,cax)
    axis xy
    colorbar
    title(sprintf('Ratio C_{lm}: h_{lm}/g_{lm}, should be %g',V(index)))
    xlabel('m')
    ylabel('l')  
    
    axes(ah(3))
    cax=V(index)*[0.75 1.25];
    imagesc([-L,L],[0,L],(ratio_blm+ratio_clm)/2,cax)
    axis xy
    colorbar
    title(sprintf(...
      'Ratio average C_{lm},B_{lm}: h_{lm}/g_{lm}, should be %g',V(index)))
    xlabel('m')
    ylabel('l')  
    
    nans=find(isnan(ratio_blm)); 
    ratio_blm(nans)=0;
    nans=find(isnan(ratio_clm));
    ratio_clm(nans)=0;
    avgratio_blm=sum(sum(abs(ratio_blm)))/(L+1)/(L+1);
    avgratio_clm=sum(sum(abs(ratio_clm)))/(L+1)/(L+1);
    disp(sprintf('Average ratio Blm = %g, Clm = %g, total = %g',...
        avgratio_blm,avgratio_clm,(avgratio_blm+avgratio_clm)/2));   
    figdisp('spectral_ratio',sprintf('%s_%s_%i_%i',comp,dom,L,index))
    
end

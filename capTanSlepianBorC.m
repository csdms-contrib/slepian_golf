function [GB,V]=capTanSlepianBorC(L,TH)
  % [GB,V]=capTanSlepianBorC(L,TH)
  %
  % Calculates the Slepian transformation matrix if we are only interested
  % in the concentration of either the Blm or Clm
  %
  % INPUT:
  %
  % L    maximum spherical-harmonic degree
  % TH   spherical cap opening angle
  %
  % OUTPUT:
  %
  % GB   Slepian transformation matrix
  % V    concentration factors
  %
  % Last modified by plattner-at-alumni.ethz.ch, 11/19/2018
  


  
if ~ischar(L)
  
  % Fromvectanglmalpha.m
  mvec=0:L;
  sizesBC=max(L+1-max(mvec,1), zeros(size(mvec)) );
  siztot=sizesBC;%2*sizesBC;
  
  deM=addmout(L);
  
  alpha=cumsum([1 siztot(1) gamini(siztot(2:end),2) ]);
  
  nelems=(L+1)^3 - L*(L+1)^2 +L*(L+1)*(2*L+1)/6;
  %GB=sparse((L+1)^2,2*(L+1)^2 - 2);
  GB=sparse([],[],[],(L+1)^2,(L+1)^2 - 1, nelems);
  % Treat Blm,Clm like Plm and remove the zero part later

  V=nan(1,(L+1)^2-1);
  
  parfor mm=1:L+1
    m=mm-1;
    [~,~,Bm]=kerneltancapm(TH,L,m);
    [Cm,Vm]=eig(Bm);
    [Vm,isrtm]=sort(sum(Vm,1),'descend');
    Cm=Cm(:,isrtm);
    Vppos{mm}=Vm;
    Vpneg{mm}=Vm;
    CBpos{mm}=Cm;
    CBneg{mm}=Cm;
  end

  CBpos{1} = [zeros(1,size(CBpos{1},2)); CBpos{1}];
  CBneg{1} = [zeros(1,size(CBpos{1},2)); CBpos{1}];

  for m=0:L
    if m>0
      GB(deM==-m,alpha(2*m):alpha(2*m+1)-1)=CBneg{m+1};
      V(alpha(2*m):alpha(2*m+1)-1)=Vpneg{m+1};
    end
    %keyboard
    GB(deM==m,alpha(2*m+1):alpha(2*m+2)-1)=CBpos{m+1};
    V(alpha(2*m+1):alpha(2*m+2)-1)=Vppos{m+1};
  end

  GB=GB(2:end,:);

  % Now sort
  [V,isrt]=sort(V,'descend');
  GB=GB(:,isrt);
  

  elseif strcmp(L,'demo')
    index=3;
    L=20;
    TH=20;
    [GB,V]=capTanSlepianBorC(L,TH);
    % Assume it's Blm or Clm
    blmcosi=fcoef2flmcosi(GB(:,index),1);
    clmcosi=blmcosi;
    %clmcosi(:,3:4)=zeros(size(clmcosi(:,3:4)));
    blmcosi(:,3:4)=zeros(size(blmcosi(:,3:4)));

    % Eval and plot
    [r,lon,lat]=blmclm2xyz(blmcosi,clmcosi,1);
    caxis([-1,1]*max(abs(caxis)))
    kelicol(1)
    subplot(1,2,1)
    title('colatitudinal component')
    plotplm(r(:,:,2),lon*pi/180,lat*pi/180,2);
    view(90,90)
    subplot(1,2,2)
    title('longitudinal component')
    plotplm(r(:,:,1),lon*pi/180,lat*pi/180,2);
    view(90,90)
    caxis([-1,1]*max(abs(caxis)))
    kelicol(1)
end

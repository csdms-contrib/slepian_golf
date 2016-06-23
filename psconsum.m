function psconsum(L)
% PSCONSUM(L)
%
% Plots the eigenvalue-weighted sum of all vector Slepian fields localized 
% to all of the continents for the radial and the tangential components.
% As in: Simons, Dahlen and Wieczorek (SIAM, 2006), Figure 6.4
%
% INPUT:
%
% L       Spherical harmonic degree of the bandwidth [default: 18]
%
% Last modified by plattner-at-alumi.ethz.ch, 3/1/2012
%
% See also PSALLCONS

doms={'africa', 'eurasia', 'namerica', 'australia', 'greenland', ...
      'samerica','antarctica'};
sa=0;
for index=1:length(doms)
  sa=sa+spharea(doms{index});
end
defval('L',18);

% Calculate the Shannon number
N=( (L+1)^2 + 2*((L+1)^2-1) )*sa; % Size of radial + tangential space

% Specify the truncation levels in the sum
ens=ceil([N/4 N/2 N 3*(L+1)^2-2]);

legsi{1}=sprintf('1%s N/4','\rightarrow');
legsi{2}=sprintf('1%s N/2','\rightarrow');
legsi{3}=sprintf('1%s N','\rightarrow');
legsi{4}=sprintf('1%s 3(L+1)^2-2','\rightarrow');

clf
[ah,ha]=krijetem(subnum(2,2));


for index=1:length(ens)
  axes(ah(index))
  fnpl=sprintf('%s/PSALLCONS-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'PSALLCONS'),L,ens(index));
  if exist(fnpl,'file')==2
    load(fnpl)
  else
    F=psallcons(doms,L,0);
  end
  % The maximum value is smaller or equal to N/A
  NA=3*(L+1)^2-2;
  F(F<NA/100)=NaN;
  imagefnan([0 90],[360 -90],F,'kelicol',[-NA NA])
  set(gca,'ytick',[-90:45:90])
  if index==3
    set(gca,'xtick',[0:90:270])
  else
    set(gca,'xtick',[0:90:360])
  end
  deggies(gca)
  [jk,a]=plotcont; set(a,'linew',1)
  axis image  
  [bh(index),th(index)]=boxtex('ll',ah(index),legsi{index},14);
end

% Cosmetics
longticks(ah,3/2)
nolabels(ha(3:4),2)
nolabels(ah(1:2),1)
set(ah,'Camerav',6.5)
serre(ha(1:2),1.25,'down')
serre(ha(3:4),1.25,'down')

% TRICK TO GET A COLOR BAR
caxcon=[0 1];
caxoc=[-1 0];
h=axes;
[cb,xcb]=addcb('hor',caxcon,caxoc,'kelicol(1)');
delete(h)
shrink(cb,2,2)
movev(cb,-0.075)
set(xcb,'String','cumulative energy')
axes(cb)
xlim([0 1])

fs=15;
set(cb,'xtick',[0 1],'xtickl',{'0' 'N/A'},'FontS',fs)
set(xcb,'FontS',fs)
set(ah,'FontS',fs-1)

movev([ah cb],.1)

fig2print(gcf,'landscape')
figdisp

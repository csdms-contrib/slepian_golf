function psvals(all)
% PSVALS(all)
%
% Eigenvalue structure for various concentration regions.
% Mixed-order ranking. Tangential. Show one symbol per eigenvalue or one
% symbol for every double-couple of eigenvalues
%
% INPUT:
%
% all   Show all eigenvalues (1) or only one per couple (0)
% 
% Last modified by plattner-at-alumni.ethz.ch, 2/29/2012
%
% See also PSELM

defval('all',1)
TH=[10 20 30 40];
L=18;

clf
[ah,ha,H]=krijetem(subnum(2,2));

yls=[-0.1 1.1];

symbs={'o','x','s','+','v','*','^','d','<','>','p','h',...
       'o','x','s','+','v','*','^','d'};
xmax=120;

more off

for index=1:length(TH)
  legsi{index}=sprintf('%s = %i%s ','\Theta',TH(index),str2mat(176));
  Nall(index)=(2*(L+1)^2-2)*(1-cos(TH(index)*pi/180))/2;

  [lrnk,mrnk,lval,VV,Vsum]=pselm(TH(index),L);
  
  ldubs=lval;
  all=mod(all,2)
  step=2-all;
  axes(ah(index))
  for ondi=1:step:length(ldubs)
    p(ondi,index)=plot(ondi,ldubs(ondi),symbs{mrnk(ondi)+1});
    hold on
  end
  % plot(lval,'--r')
  hold on
  plot(round([Nall(index) Nall(index)]),yls,'k:')
  plot([0 xmax],[0.5 0.5],'k:')
  plot([0 xmax],[0 0],'k:')
  plot([0 xmax],[1 1],'k:')
  set(ah(index),'xlim',[0 xmax],'ylim',yls,'xgrid','off','ygrid','off',...
		'xticklabel',[1 10:20:xmax],'xtick',[1 10:20:xmax],...
		'ytick',[0:0.25:1])
  try 
    [bh(index),th(index)]=boxtex('ur',ah(index),legsi{index},...
          12,[],[],[],0.8);
  end
  drawnow
end

axes(ah(4))
fact=1.6;
fb=fillbox([2/fact fact*18 0.88 -0.05],'w');
% m=0 is not pm
ondi=1;
ypo=0+0.075*(ondi-1);
pl(ondi,1)=plot(4,ypo,symbs{ondi});
hold on
tl(ondi,1)=text(7,ypo,sprintf(' m = %s %i','  ',ondi-1),'FontSize',8);

for ondi=2:12
  ypo=0+0.075*(ondi-1);
  pl(ondi,1)=plot(4,ypo,symbs{ondi});
  hold on
  tl(ondi,1)=text(7,ypo,sprintf(' m = %s %i','\pm',ondi-1),'FontSize',8);
end

% Now make the plot beautiful
try 
    movev(th,-.00)
end
longticks(ah)
set([p(~~p(:)) ; pl(~~pl(:))],'MarkerSize',4,'MarkerFaceColor',grey,'MarkerEdgeColor','k')
axes(ha(1))
al(1)=ylabel('eigenvalue \lambda');
axes(ha(2))
al(2)=ylabel('eigenvalue \lambda');
axes(ha(2))
xl(1)=xlabel('rank');
axes(ha(4))
xl(2)=xlabel('rank');

nolabels(ha(3:4),2)
nolabels(ah(1:2),1)

serre(H',1/2,'down')
serre(H,1/2,'across')

for ind=1:4
  xx(ind)=xtraxis(ah(ind),round(Nall(ind)),...
		  {sprintf('N = %i',round(Nall(ind)))});
end
longticks(xx)

set([xl al],'FontSize',13)
set([ ah],'FontSize',12)

figdisp



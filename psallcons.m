function [F,Vfull,K]=psallcons(dom,L,mkmov)
%[F,V,K,Vc,Kc]=psallcons(dom,L,mkmov)
% [F,V,K,Vc,Kc]=PSALLCONS(dom,L,mkmov)
%
% Similar to Simons, Dahlen and Wieczorek (2006): Localization to the 
% continents but for vector fields on the sphere, both, the radial and the
% tangential component
%
% INPUT:
%
% dom     Cell with domains (strings) you want included [default: all]
% L       Maximal spherical harmonic bandwidth [default: 18]
% mkmov   1 Make a movie of this
%         0 Don't [default]
%
% OUTPUT:
%
% F       The eigenvalue-weighted sum of the energies
% Vfull   Ordered list of eigenvalues, component (radial (0) or tangential
%         (1)) and rank within component 
% K       Localization kernel
%
% Last modified by plattner-at-alumi.ethz.ch, 3/1/2012
%
% See also PSCONSUM, KERNELBP, KERNELCP

defval('dom',...
       {'africa','eurasia','namerica','australia','greenland', ...
	'samerica','antarctica'});
defval('L',18);
defval('mkmov',0);
degres=0.5;

if ~iscell(dom)
  dom={dom};
end

%% Start with the radial part
rotb=0;
% Initialize the kernel matrix
area=0;
Krad=repmat(0,(L+1)^2,(L+1)^2);
for index=1:length(dom)
  area=area+spharea(dom{index});  
  disp(sprintf('Working on %s',dom{index}))
  if strcmp(dom{index},'antarctica')
      rotb=1;
  else
      rotb=0;
  end
  [Klmlmp,XY]=kernelcp(L,dom{index},[],[],rotb);
  Krad=Krad+Klmlmp;
  clear Klmlmp
end

% Stolen from LOCALIZATION
[Crad,Vrad]=eig(Krad);
[Vrad,isrt]=sort(sum(real(Vrad),1));
Vrad=fliplr(Vrad);
Crad=Crad(:,fliplr(isrt));

% Sticks the cosine/sine coefficients back
% into the right place of LMCOSI
[dems,dels,mz,lmc,mzin]=addmon(L);
for index=1:size(Crad,2)
  CCrad{index}=reshape(insert(Crad(:,index),0,mzin),2,length(dems))';
end

% Eigenvalues and eigenfunctions
Vrad=Vrad(:);
Crad=CCrad;



%% Now the tangential part
% Initialize the kernel matrix
Ktan=repmat(0,2*(L+1)^2-2,2*(L+1)^2-2);
for index=1:length(dom)
  disp(sprintf('Working on %s',dom{index}))
  if strcmp(dom{index},'antarctica')
      rotb=1;
  else
      rotb=0;
  end
  [Klmlmp,XY]=kernelbp(L,dom{index},[],[],rotb);
  Ktan=Ktan+Klmlmp;
  clear Klmlmp
end

% Stolen from VECTORSLEPIAN
[Ctan,Vtan]=eig(Ktan);
[Vtan,isrt]=sort(sum(real(Vtan),1));
Vtan=fliplr(Vtan);
Ctan=Ctan(:,fliplr(isrt));

CBlm=Ctan(1:size(Ctan,1)/2,:);
CClm=Ctan((size(Ctan,1)/2+1):end,:);  

CBlm=[NaN(1,size(Ctan,2));CBlm]; 
CClm=[NaN(1,size(Ctan,2));CClm];  

% Sticks the cosine/sine coefficients back
% into the right place of LMCOSI
% Make sure that the l=0 component is left out (that's why it is more 
% complicated) 
[demsz,delsz,mz,lmc,mzin]=addmon(L);
dems=demsz(2:end);
dels=delsz(2:end);
for index=1:size(Ctan,2)
  CCB=reshape(insert(CBlm(:,index),0,mzin),2,length(demsz))';
  CCC=reshape(insert(CClm(:,index),0,mzin),2,length(demsz))';     
  CCBlm{index}=CCB(2:end,:);  
  CCClm{index}=CCC(2:end,:);     
end

% Eigenvalues and eigenfunctions
Vtan=Vtan(:);
CBlm=CCBlm;
CClm=CCClm;


% Make a full mixed ranking of the radial and tangential Slepians
% First col: the actual eigenvalues
% Second col: if they are radial (0) or tangential (1)
% Third col: The rank inside radial or tangential 
Vfull=[Vrad zeros(size(Vrad)) (1:length(Vrad))';Vtan ones(size(Vtan)) (1:length(Vtan))'];

Vfull=sortrows(Vfull);
Vfull=flipud(Vfull);



%% Now wrap everything up
% Provide the progressive spatial sum
F=0;

if mkmov==1
  % MAKE SURE THE FIGURE FILE IS UNINTERRUPTED VIEW
  mov=avifile('psmovie.avi');
end
for index=1:3*(L+1)^2-2
  typ=Vfull(index,2);
  innerindex=Vfull(index,3);
  fnpl=sprintf('%s/PSALLCONS-%i-%i.mat',...
	       fullfile(getenv('IFILES'),'PSALLCONS'),L,index);
  if exist(fnpl,'file')==2
    FF=load(fnpl);
    F=FF.F;
    % This doesn't work here
    NA=3*(L+1)^2-2;
    imagefnan([0 90],[360 -90],F,[],[-NA NA],[],[],100);
    plotcont
    title(sprintf('%3.3i',index)) 
    axis tight
    if mkmov==1
      % Make and save movie
      mov=addframe(mov,getframe(gcf));
    end
  else    
      if(typ==0)
        F=F+(plm2xyz([delsz demsz Crad{innerindex}],degres)).^2*Vrad(innerindex);
        save(fnpl,'F')
      else 
        data=blmclm2xyz([dels dems CBlm{innerindex}], ...
            [dels dems CClm{innerindex}],degres);
        F=F+(data(:,:,1).^2+data(:,:,2).^2)*Vtan(innerindex);  
        save(fnpl,'F')
      end
  end
end
N=sum(Vrad)+sum(Vtan);
disp(sprintf('Max function value is %f, N/A is %f, 3*(L+1)^2)-2 = %f ',...
    max(max(F)),N/area,3*(L+1)^2-2));
if mkmov==1
  mov=close(mov);
end


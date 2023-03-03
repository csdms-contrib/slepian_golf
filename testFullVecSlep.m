function testFullVecSlep(index)
% Use this to test if glmalphaFullVec does everything correctly
% Last modified by plattner-at-alumni.ethz.ch, 2023/3/2

figshift = 0;%3;
Lmax = 110;
%dom = 'namerica';
dom = 10;
rotcoord = [35,90];
%rotcoord = [];
J=100;
onorout = 0; 


deg = 1;

[G,V] = glmalphaFullVec(Lmax,dom,1,rotcoord,J,onorout);

    % Evaluate Elm
    coefE = G(1:(Lmax+1)^2 ,index);
    [valE,lon,lat]=elm2xyz(coef2lmcosi(coefE,onorout), deg);

    % Evaluate Flm
    coefF = G((Lmax+1)^2+1:2*(Lmax+1)^2-1,index);
    [valF,lon,lat]=flm2xyz(fcoef2flmcosi(coefF,onorout), deg);

    % Evaluate Clm
    coefC = G(2*(Lmax+1)^2:end,index); 
    coefB = zeros(size(coefC));
    [valC,lon,lat] = blmclm2xyz(fcoef2flmcosi(coefB,onorout), fcoef2flmcosi(coefC,onorout),deg);

    % Add them up
    radialEFC = valE{1} + valF{1};
    longitudinalEFC = valE{3} + valF{3} + valC(:,:,1);
    colatitudinalEFC = valE{2} + valF{2} + valC(:,:,2);



    % Plot the components. This should be the same as for the PBC basis
    figure(1+figshift)

    subplot(3,1,1)
    imagesc(lon,lat,radialEFC)
    title('radial component EFC')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
    colorbar

subplot(3,1,2)
    imagesc(lon,lat,longitudinalEFC)
    title('longitudinal component EFC')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
    colorbar

subplot(3,1,3)
    imagesc(lon,lat,colatitudinalEFC)
    title('colatitudinal component')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
colorbar


clear G
clear V


    [G,V] = glmalphaFullVec(Lmax,dom,0,rotcoord,J,onorout);


    

    % Plot radial component
    figure(2+figshift)

    coefP = G(1:(Lmax+1)^2,index);
    [radial,lon,lat] = plm2xyz(coef2lmcosi(coefP,onorout),deg);
    subplot(3,1,1)
    imagesc(lon,lat,radial)
    title('radial component')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
    colorbar

    % plot colatitudinal components
    coefB = G((Lmax+1)^2+1:2*(Lmax+1)^2-1,index);
    coefC = G(2*(Lmax+1)^2:end,index);
    [tangential,lon,lat] = blmclm2xyz(fcoef2flmcosi(coefB,onorout), fcoef2flmcosi(coefC,onorout),deg);
    subplot(3,1,2)
    imagesc(lon,lat,tangential(:,:,1))
    title('longitudinal component')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
    colorbar

    subplot(3,1,3)
    imagesc(lon,lat,tangential(:,:,2))
    title('colatitudinal component')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
    colorbar



figure(3+figshift)
     subplot(3,1,1)
    imagesc(lon,lat,radialEFC - radial)
    title('radial component difference')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
    colorbar

subplot(3,1,2)
    imagesc(lon,lat,longitudinalEFC - tangential(:,:,1))
    title('longitudinal component difference')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
    colorbar

subplot(3,1,3)
    imagesc(lon,lat,colatitudinalEFC - tangential(:,:,2))
    title('colatitudinal component difference')
    axis xy
    caxis([-1,1]*max(abs(caxis)));
    kelicol(1)
colorbar
function [dlon,dcola,dr]=dcart2dsph(lon,cola,dx,dy,dz)
% [dlon,dcola,dr]=dcart2dsph(lon,cola,dx,dy,dz)
%
% Calculates the vectorial components in cartesian coordinates into 
% vectorial components in spherical coordinates. Orthogonality is achieved
% by using position on the unit sphere. The transformation matrix is
% automatically orthogonal.
%
% INPUT:
%
% lon       longitudinal position (0<=lon<2*pi)
% cola      colatitudinal position (0<=lat<=pi)
% dx,dy,dz  components of the vector
%
% OUTPUT
% dlon      longitudinal component of the vector 
% dcola     colatitudinal component of the vector 
% dr        radial component of the vector
%
% Last modified by plattner-at-alumni.ethz.ch, 1/29/2014
% fixed a typo in help page on 3/7/2017

lonsize=size(lon);
lon=lon(:);
cola=cola(:);
dx=dx(:);
dy=dy(:);
dz=dz(:);

dlon=nan(size(lon));
dcola=nan(size(lon));
dr=nan(size(lon));

for i=1:length(lon)
    % See e.g. http://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
     D=[sin(cola(i))*cos(lon(i)) sin(cola(i))*sin(lon(i))  cos(cola(i));
       cos(cola(i))*cos(lon(i)) cos(cola(i))*sin(lon(i)) -sin(cola(i));
       -sin(lon(i))             cos(lon(i))                0           ];
    
    dsph=D*[dx(i);dy(i);dz(i)];
    dr(i)=dsph(1);
    dcola(i)=dsph(2);
    dlon(i)=dsph(3);
end

dlon=reshape(dlon,lonsize(1),lonsize(2));
dcola=reshape(dcola,lonsize(1),lonsize(2));
dr=reshape(dr,lonsize(1),lonsize(2));
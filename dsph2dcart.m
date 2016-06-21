function [dx,dy,dz]=dsph2dcart(lon,cola,dlon,dcola,dr)
% [dx,dy,dz]=dsph2dcart(lon,cola,dlon,dcola,dr)
%
% Calculates the vectorial components in spherical coordinates into 
% vectorial components in cartesian coordinates. Orthogonality is achieved
% by using position on the unit sphere. The transformation matrix is
% automatically orthogonal.
%
% INPUT:
%
% lon    longitudinal position of the vector (0<=lon<2*pi)
% cola   latitudinal position of the vector (0<=cola<=pi)
% dlon   longitudinal component of the vector
% dcola  colatitudinal component of the vector 
% dr     radial component of the vector
% 
% OUTPUT: 
%
% dx     x-component of the vector
% dy     y-component of the vector
% dz     z-component of the vector
%
% Last modified by plattner-at-alumni.ethz.ch, 1/29/2014

lonsize=size(lon);
lon=lon(:);
cola=cola(:);
dlon=dlon(:);
dcola=dcola(:);
dr=dr(:);

dx=nan(size(lon));
dy=nan(size(lon));
dz=nan(size(lon));

for i=1:length(dlon)
    % See e.g. http://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates
    D=[sin(cola(i))*cos(lon(i)) sin(cola(i))*sin(lon(i))  cos(cola(i));
       cos(cola(i))*cos(lon(i)) cos(cola(i))*sin(lon(i)) -sin(cola(i));
       -sin(lon(i))             cos(lon(i))                0           ];
   
   dcart=D'*[dr(i);dcola(i);dlon(i)];
   dx(i)=dcart(1);
   dy(i)=dcart(2);
   dz(i)=dcart(3);
    
end

dx=reshape(dx,lonsize(1),lonsize(2));
dy=reshape(dy,lonsize(1),lonsize(2));
dz=reshape(dz,lonsize(1),lonsize(2));
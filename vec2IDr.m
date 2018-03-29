function [I,D,r]=vec2IDr(dr,dcola,dlon)
% [I,D,r]=vec2IDr(dr,dcola,dlon)
%
% Calculates inclination, declination, and length from longitudinal, 
% latitudinal, and radial component of a vector
%
% INPUT:
%
% dr       radial component (outward)
% dcola    colatitudinal component (southward)
% dlon     longitudinal component (eastward)
%
% OUTPUT:
%
% I        inclination [radians]
% D        declination [radians]
% r        length
%
% Last modified by plattner-at-alumni.ethz.ch, 3/29/2018

r = sqrt(dlon.^2 + dcola.^2 +dr.^2);
h = sqrt(dlon.^2 + dcola.^2);
D = atan2(-dcola,dlon)+pi;
I = atan(-dr./h);




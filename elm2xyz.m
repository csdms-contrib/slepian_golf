function varargout=elm2xyz(elmcosi,degrees,c11cmn,lmax)
% [r,lon,lat]=ELM2XYZ(elmcosi,degres,c11cmn,lmax)
%
% Inverse inner source vector spherical harmonic transform
%
% Compute a gradient spherical harmonic vector field (coefficents for 
% r, theta, and phi) from coefficients of the Elm given as [l m Ccos Csin] 
% (not necessarily starting from zero, but sorted) with degree resolution 
% 'degres' [default: approximate Nyquist degree].
%
% Using 4*pi-normalized real spherical harmonics
% INPUT:
%
% elmcosi   Matrix listing l,m,cosine and sine expansion coefficients for
%           Elm
% degres    Longitude/ latitude spacing, in degrees [default: Nyquist] OR
%           "lat": a column vector with latitudes [degrees]
% c11cmn    Corner nodes of lon/lat grid [default: 0 90-sqrt(eps) 360 -90]
%           OR "lon": a column vector with longitudes [degrees]
% lmax      Maximum bandwidth expanded at a time [default: 720]
%
% OUTPUT: 
%
% r         The field (matrix for a grid, vector for scattered points)
%           r{1} is the radial component; r{2} is the theta-component,
%           r{3} is the phi component.
%           First dimension of r{i} is lat, second dimension is lon.
% lon,lat   The grid (matrix) or evaluation points (vector), in degrees
%
% See also BLMCLM2XYZ, PLM2XYZ
%
% Last modified by plattner-at-alumni.ethz.ch 02/24/2015

warning('Using 4*pi-normalized real spherical harmonics')
  
defval('degrees',1)
defval('c11cmn',[])
defval('lmax',[])
tol=1e-5; % Tolerance for testing if the positions for the longitudes and 
          % latitudes are the same 

% Start with transforming elm into plm & blm.
% Generate flmcosi. Remember it has no zero degree
flmcosi=[elmcosi(2:end,1:2) zeros(size(elmcosi(2:end,3:4)))];
[plmcosi,blmcosi]=EFC2PBC(elmcosi,flmcosi);
% Then apply blmclm2xyz
% First make a clm. 
clmcosi=[blmcosi(:,1:2) zeros(size(blmcosi(:,3:4)))];
[rtan,lont,latt]=blmclm2xyz(blmcosi,clmcosi,degrees,c11cmn,lmax);
% Then apply plm2xyz
[rrad,lonr,latr]=plm2xyz(plmcosi,degrees,c11cmn,lmax);
% then put the plm result into the radial component and the blmclm result
% into the tangential component
r{1}=rrad;
r{2}=rtan(:,:,2);
r{3}=rtan(:,:,1);
% Test if something is wrong with the lon and lat positions
if(norm(lont-lonr)>tol || norm(latt-latr)>tol)
    error('something is wring with the positions')
end
vars={r,lont,latt};
varargout=vars(1:nargout); 


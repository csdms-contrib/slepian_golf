function varargout=preparetanplot(data)
% [dat,dmax]=PREPARETANPLOT(data)
%
% Returns the intensity (if data is a struct, then the intensity of 
% data{1})
%
% INPUT:
%
% data  data(:,:,1) is a 2d array containing the values of the vector 
%       field in longitudinal direction,(phi) 
%       data(:,:,2) is a 2d array containing the values of the vector 
%       field in latitudinal direction (theta)
%       OR:
%       data is a struct, where data{1} as mentioned before
% 
% OUTPUT:
%
% dat   2d array containing the intensity values of the vectors (the 
%       negative to make the plotting colors in kelicol red). First 
%       dimension is lat, second is lon)
% dmax  Maximum intensity value
%
% See also KELICOL
%
% Last modified by plattner-at-alumni.ethz.ch, 02/27/2012

if iscell(data)
     absdata=sqrt(data{1}(:,:,1).^2+data{1}(:,:,2).^2);
else
     absdata=sqrt(data(:,:,1).^2+data(:,:,2).^2);
end

dmax=max(max(absdata)); thresh=100;  
dat=-absdata;
dat(abs(dat)<dmax/thresh)=0;

varns={dat,dmax};
varargout=varns(1:nargout);
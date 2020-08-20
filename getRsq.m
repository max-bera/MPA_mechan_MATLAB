function Rsq=getRsq(data)

% -------------------------------------------------------------------------
% getRsq.m - calculates R-squared value
%
% Rsq=getRsq(data)
%
% input:
% data      - n-by-2 matrix with    data(:,1) - observed data
%                                   data(:,2) - expected data
%
% output:
% Rsq        - Resulting R-2 value
%
% subfunctions:
% none
% -------------------------------------------------------------------------
if nargin<1, help getRsq; end

SStot=sum((data(:,1)-mean(data(:,1))).^2);
SSfit=sum((data(:,1)-data(:,2)).^2);
Rsq=1-SSfit/SStot;
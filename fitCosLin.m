function [omega,phi,ampl,y0,LinFit]=fitCosLin(x,y,per,forceper,initLin,pl)

% -------------------------------------------------------------------------
% fitCosLin.m - fit cosine with superposed linear function 
% 
% [omega,phi,ampl,y0,LinFit]=fitCosLin(x,y,per,forceper,initLin,pl)
%
%
% input:
% x         - x-values for sinusoidal function
% y         - y-values for sinusoidal function
% per       - approximate period (optional)
% forceper  - toggle forced periodicity (1=periodicity is per, 0=fits per
%               (default))
% initLin   - initial value for linear fit ([offset slope])
% pl        - toggle plot
%
% output:
% omega     - angular frequency (omega = 2*pi / T)
% phi       - phase shift 
% ampl      - amplitude 
% y0        - zero offset
% LinFit    - linear fit result
% 
% Corresponding to: y = linFit(1) + linFit(2) * x + ampl * cos(omega * x + phi)
%
% subfunctions
% all cosine-difference functions are at the bottom of this m-file
% -------------------------------------------------------------------------
if nargin<1, help fitCos; return; end
if nargin<6, pl=0; end

if size(x,1)==1; x=x'; end
if size(y,1)==1; y=y'; end

if nargin<3 || isempty(per)
    % obtain periodicity using fft
    incpp   = mean((x(2:end)-x(1))./[1:length(x)-1]');  % increase per datapoint
    Fs      = 1/incpp; % sampling frequency
    [Y,f]   = periodogram(y,[],[],Fs);
    % L       = length(x);
    % NFFT    = 2^nextpow2(L); % next power of 2 from length of y
    % Y       = abs(fft(y,NFFT))/L; % fourier spectrum
    % f       = Fs/2*linspace(0,1,NFFT/2+1);
    f(1:2)=[]; Y(1:2)=[]; % cut off center peak
    fmax    = f(Y==max(Y)); % peak location
    per     = abs(1/fmax);   % approximate period
end
if nargin<4 || isempty(forceper), forceper=0; end

% get starting conditions for phase-lag fit and fit phase only (!)
y0      = (max(y)+min(y))/2;  % estimation offset
ampl0    = (max(y)-min(y))/2; % calc amplitude
% phi0    = 2*pi*mod(x(max(y)==y),per)/per;

% find the approximate phaseshift
phi0    = pi/2;
phi     = fminsearch(@(phi) cosDiff_phi([x,y],ampl0,per,phi),phi0);

% perform final fit on both amplitude, phase and (optional) periodicity
vars0   = [ampl0,phi,per];
if forceper
    vars0   = [vars0,initLin];
    options = optimset('MaxFunEvals',5e3)%,'PlotFcns',@optimplotfval)
	vars    = fminsearch(@(vars) cosDiff_forceperLinFit([x,y],vars,per),vars0,options);
    
else
    vars    = fminsearch(@(vars) cosDiff([x,y],vars),vars0);
    vars0   = [vars,initLin];
    vars    = fminsearch(@(vars) cosDiff_forceperLinFit([x,y],vars,vars0(3)),vars0);
    per     = vars(3); 
end
ampl    = vars(1);
phi     = vars(2);
omega   = 2*pi/per;
LinFit  = vars(4:5);

if pl
	hF=findobj('type','figure','name','fitCosLin');
    if isempty(hF)
        hF=figure;
        set(hF, 'name', 'fitCosLin');
    end
    figure(hF);clf
    set(hF,'Color',[1 1 1]);
    
    plot(x,y,'.');
    hold on
    xfit=linspace(min(x),max(x),1e3);
	yfit=ampl.*cos(omega.*xfit+phi)+LinFit(1)+LinFit(2).*xfit;
    plot(xfit,yfit,'-');
end
end

% additional functions for fminsearch use:
function diff=cosDiff_phi(data,V,per0,phi)

yfit    = V.*cos((2*pi/per0).*data(:,1)+phi); 
diff    = sum((yfit-data(:,2)).^2);
end
function diff=cosDiff_forceper(data,vars,per)

yfit    = vars(1).*cos((2*pi/per).*data(:,1)+vars(2)); 
diff    = sum((yfit-data(:,2)).^2);
end
function diff=cosDiff_forceperLinFit(data,vars,per)

yfit    = vars(1).*cos((2*pi/per).*data(:,1)+vars(2))+vars(4)+vars(5).*data(:,1); 
diff    = sum((yfit-data(:,2)).^2);
end
function diff=cosDiff(data,vars)

yfit    = vars(1).*cos((2*pi/vars(3)).*data(:,1)+vars(2)); 
diff    = sum((yfit-data(:,2)).^2);
end
function diff=cosDiff_LinFit(data,vars)

yfit    = vars(1).*cos((2*pi/vars(3)).*data(:,1)+vars(2))+vars(4)+vars(5).*data(:,1); 
diff    = sum((yfit-data(:,2)).^2);
end
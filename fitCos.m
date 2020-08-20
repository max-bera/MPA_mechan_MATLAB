function [omega,phi,ampl,y0,LinFit]=fitCos(x,y,per,forceper,pl)

% -------------------------------------------------------------------------
% fitCos.m - fit cosine function to sinusoidal data
% 
% [omega,phi,ampl,y0]=fitCos(x,y,per,forceper,pl)
%
% NB: this function fits a sinusoidal signal withthese steps:    
%           (1) fft for period estimation (unless per is given as input)
%           (2) phase fit with ampl and offset guess from min and max
%           (3) fit all together period, phase and amplitude
%           (4) fit offset with previously determined period phase and ampl
%
% input:
% x         - x-values for sinusoidal function
% y         - y-values for sinusoidal function
% per       - approximate period (optional, but highly recommended)
% forceper  - toggle forced periodicity (1=periodicity must be per, 
%               0=fits period - default)
% pl        - toggle plot (1=on, 0=off - def)
% 
% output:
% omega     - angular frequency (omega = 2*pi / period)
% phi       - phase shift (radians)
% ampl      - amplitude 
% y0        - offset
% 
% Corresponding to: y = y0 + ampl * cos(omega * x + phi)
%
% subfunctions
% all cosine-difference functions are at the bottom of this m-file
% -------------------------------------------------------------------------
if nargin<1, help fitCos; return; end
if nargin<5, pl=0; end

if size(x,1)==1; x=x'; end
if size(y,1)==1; y=y'; end

if nargin<3 || isempty(per) % if no period is given in, then use fft to guess it
    % obtain periodicity using fft
    incpp   = mean((x(2:end)-x(1))./[1:length(x)-1]');  % increase per datapoint
    Fs      = 1/incpp; % sampling frequency
    [Y,f]   = periodogram(y,[],[],Fs);
    % L       = length(x);
    % NFFT    = 2^nextpow2(L); % next power of 2 from length of y
    % Y       = abs(fft(y,NFFT))/L; % fourier spectrum
    % f       = Fs/2*linspace(0,1,NFFT/2+1);
    f(1:3)=[]; Y(1:3)=[]; % cut off center peak
    fmax    = f(Y==max(Y)); % peak location
    per     = abs(1/fmax);   % approximate period
%     figure;plot(f,Y)
end
if nargin<4 || isempty(forceper), forceper=0; end

% get starting conditions for phase-lag fit and fit phase only (!)
y0      = (max(y)+min(y))/2;  
ampl0   = (max(y)-min(y))/2; % calc amplitude
% phi0    = 2*pi*mod(x(max(y)==y),per)/per;

% find the approximate phaseshift
phi0    = pi/2;
phi     = fminsearch(@(phi) cosDiff_phi([x,y],y0,ampl0,per,phi),phi0);

% perform fit on both amplitude, phase and (optional) periodicity
if forceper
    % in this option the period given by per is forced on the data, only
    % amplitude and phase are fitted
    vars0   = [ampl0,phi];
    vars    = fminsearch(@(vars) cosDiff_forceper([x,y],vars,per,y0),vars0);
    vars(3) = per;
else
    % in this option the period is fitted with starting value per next to
    % the variables amplitude and phase
    vars0   = [ampl0,phi,per];
    vars    = fminsearch(@(vars) cosDiff([x,y],vars,y0),vars0);
    per     = vars(3);
end
% finally the offset is fitted to the data with the previous variables
y0      = fminsearch(@(y0i) cosDiff_y0([x,y],vars,y0i),y0);

% generate output
ampl    = vars(1);
phi     = vars(2);
omega   = 2*pi/per;

if pl
	hF=findobj('type','figure','name','fitCos');
    if isempty(hF)
        hF=figure;
        set(hF, 'name', 'fitCos');
    end
    figure(hF);clf
    set(hF,'Color',[1 1 1]);
    
    plot(x,y,'.');
    hold on
    xfit=linspace(min(x),max(x),1e3);
    yfit=y0+ampl.*cos(omega.*xfit+phi);
    plot(xfit,yfit,'-');
end
end

% additional functions for fminsearch use:
function diff=cosDiff_phi(data,y0,V,per0,phi)

yfit    = y0+V.*cos((2*pi/per0).*data(:,1)+phi); 
diff    = sum((yfit-data(:,2)).^2);
end
function diff=cosDiff_forceper(data,vars,per,y0)

yfit    = y0+vars(1).*cos((2*pi/per).*data(:,1)+vars(2)); 
diff    = sum((yfit-data(:,2)).^2);
end
function diff=cosDiff(data,vars,y0)

yfit    = y0+vars(1).*cos((2*pi/vars(3)).*data(:,1)+vars(2)); 
diff    = sum((yfit-data(:,2)).^2);
end
function diff=cosDiff_y0(data,vars,y0)

yfit    = y0+vars(1).*cos((2*pi/vars(3)).*data(:,1)+vars(2)); 
diff    = sum((yfit-data(:,2)).^2);
end

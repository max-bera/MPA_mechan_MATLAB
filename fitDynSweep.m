function [freq,E1,E2,phi,omd,phid,ampld,y0d,omf,phif,amplf,y0f]=...
            fitDynSweep(t,d,f,R,R_c,pl,ni,freq_in,fName)
%%
% fitDynSweep.m - fit one oscillation sweep and obtain results
%
%INPUT:
%t         - time (s)
%d         - aspirated length (m)
%f         - aspiration pressure (Pa)
%R         - radius of the pipette (m)
%R_c       - radius of the sample (m)
%pl        - toggle plot (1=on, 0=off)
%nu        - poisson ratio (def 0.5)
%freq_in   - frequency as starting point for sine-fit (def is determined
%            by fft of load-data input)
%fName     - filename to save results to (default not saved)
%
%output:
%freq      - frequency (Hz); fitted on pressure-sweep and imposed in
%            aspirated length sweep 
%E1        - Storage Young's Modulus E' (Pa)
%E2        - Loss Young's Modulus E'' (Pa)
%phi       - Loss tangent (degrees)
%omd       - angular frequency aspirated length oscillation
%phid      - phase fit aspirated lengthoscillation
%ampld     - amplitude fit aspirated length oscillation
%y0d       - offset fit aspirated length oscillation
%omf       - angular frequency pressure oscillation
%phif      - phase fit pressure oscillation
%amplf     - amplitude fit pressure oscillation
%y0f       - offset fit pressure oscillation
%%
if nargin<1; help fitDynSweep; return; end
if nargin<6, pl=0; end
if nargin<7, ni = 0.5; end
if nargin<8, freq_in=[]; end
if nargin<9, fName=[]; end
%%
minRsq = 0;
beta1 = 2.0142;
beta3 = 2.1187;

% offset time-data
t=t-t(1);

% fit with periodic functions
if isempty(freq_in)     
    initLin=[-0.00005, 0.00005];
    per=1/freq_in;
	[omf,phif,amplf,y0f,lin_Ffit]=fitCosLin(t,f,per,1,initLin,0);
else
    per=1/freq_in;
    %with correction
    initLin=[-0.05, 0.05];
    [omf,phif,amplf,y0f,lin_Ffit]=fitCosLin(t,f,per,1,initLin,0);
end
ffit = lin_Ffit(1) + amplf .* cos(omf .* t + phif) + lin_Ffit(2).*t;
per_force=2*pi/omf; 
% do a fit with the periodicity from the force curve imposed on the indentation
%[omd,phid,ampld,y0d]=fitCos(t,d,per_force,1);
%dfit = y0d + ampld .* cos(omd .* t + phid);

% OPTIONAL add a linear fit to the cosine fit
initLin=[-0.0000005, 0.0000005];
[omd,phid,ampld,y0d,linFit]=fitCosLin(t,d,per_force,1,initLin,0);
dfit = ampld .* cos(omd .* t + phid) +linFit(2).*t+linFit(1);

%output
phi = phif-phid; % phase lag between pressure and aspirated length
if phi<0 % if the fits just happened to switch sign, it can be 2pi further
    phi=phi+2*pi;
end
phi = mod(phi,2*pi); % the phase shift can never be more than 2pi
freq= omf/(2*pi); % frequency of indentation

%calculate storage and loss modulus [Pa] and phi [deg]
E1= 3*R/beta1*amplf/(1-(R/R_c)^beta3)/ampld*cos(phi); % E'  [Pa]
E2= 3*R/beta1*amplf/(1-(R/R_c)^beta3)/ampld*sin(phi); % E'' [Pa]
phi = phi*180/pi;

%calculate the error to the fit of the indentation depth
RSq=getRsq([d,dfit]);
 if RSq<minRsq
     E1=NaN; E2=NaN;
     phi=NaN; a=NaN;
 end

%plot results
if pl
    hF=findobj('type','figure','name','dynamicFit');
    if isempty(hF)
        hF=figure;
        set(hF, 'name', 'dynamicFit');
    end
    figure(hF);clf
    hold on
    box on
    plot(t,(ffit-lin_Ffit(1)),'--r');
    plot(t,(f-lin_Ffit(1)),'-ro','MarkerIndices',1:50:length(f));
    ylabel('\DeltaP (Pa)','Color','r');
    yyaxis right    
    plot(t,(dfit-linFit(1))*1e9,'--b')
    plot(t,(d-linFit(1))*1e9,'-b*','MarkerIndices',1:50:length(f))
    ylabel('L_p (nm)','Color','b');
    xlabel('Time (s)');
    ax=gca;
    ax.YAxis(1).Color='red';
    ax.YAxis(2).Color='blue';
    if ~isempty(fName)
        print(gcf,'-dsvg',[fName,'f',num2str(1000*freq_in),'mHz.svg']);
    end
end
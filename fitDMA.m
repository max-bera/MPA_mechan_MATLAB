function [Es,El,tand] = fitDMA(time,Lp,dP,pip_rad,Rc,freqs,ncycles,dt,...
                        rampt,creept,waitt,rate,name,plotres)
%%
%fitDMA.m takes Lp/dP and the experimental protocol as an input, segments
%the data in the single frequency portions, and calculates E'/E'' and tand
%by using an extend version of the linearized ZHou model.
%
%!!! The model currently assumens incompressibility (i.e. ni = 0.5) !!!
%
%INPUT:
%time [s]      *    time
%Lp [m]        *    array of aspirated lengths
%dP [Pa]       *    array of corresponding pressure points
%pip_rad [m]   *    pipette radius
%Rc [m]        *    specimen radius 
%freqs [Hz]    *    array of testing frequencies
%ncycles       *    # of cycles per testing frequency
%dt [s]        *    wait time before preload 
%rampt [s]     *    time to reach the preload
%creept [s]    *    creep time at target preload
%waitt [s]     *    waiting period in between oscillations
%rate {Hz]     *    sampling rate
%name               sample name (for saving *.svg plots)
%plotres            boolean, plot fit results
%
%OUTPUT:
%Es [Pa]            storage modulus
%El [Pa]            loss modulus
%tand [-]           loss factor
%%
if nargin<12; help fitDMA; return; end
if nargin<13; name = []; return; end
if nargin<14; plotres = 0; return; end
%%
%calculate start/end points in the array for each frequency
start_DMA = zeros(length(freqs),1);
start_DMA(1) = ceil(dt+rampt+creept)*rate;
fin_DMA(1) = floor(start_DMA(1)+ncycles(1)/freqs(1)*rate);
for i = 2:length(freqs)
    start_DMA(i) = floor(fin_DMA(i-1)+waitt*rate);
    fin_DMA(i) = floor(start_DMA(i)+ncycles(i)/freqs(i)*rate);
end

Es = zeros(length(freqs),1);
El = zeros(length(freqs),1);
tand = zeros(length(freqs),1);

%iterate over the frequencies, calculate parameters
for i=1:length(start_DMA)
    Lp_f = Lp(start_DMA(i):fin_DMA(i));
    dP_f = dP(start_DMA(i):fin_DMA(i));
    time_f = time(start_DMA(i):fin_DMA(i));                 
    [~,E1,E2] = fitDynSweep(time_f,Lp_f,dP_f,...
                    pip_rad,Rc,plotres,0.5,freqs(i),name);
    Es(i) = E1;
    El(i) = E2;
    tand(i) = El(i)/Es(i);
end

%plot summary of the fit
if plotres == true
    res_DMAfit = figure;
    figure(res_DMAfit)
    plot(freqs,Es/1000,'o-k')
    hold on
    xlabel('frequency [Hz]')
    plot(freqs,El/1000,'*-k')
    yyaxis left
    ylabel('E''/E'''' [kPa]');
    yyaxis right
    plot(freqs,El./Es,'.-r')
    ylabel('tand\delta');
    ax=gca;
    ax.XScale = 'log';
    ax.YColor = 'r';
    grid on
    legend E'' E'''' tan\delta;
    hold off
end
end
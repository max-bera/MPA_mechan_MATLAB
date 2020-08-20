%%
%deltasens data
close all
path = 'C:\Users\Massimiliano\Desktop\DeltaSens signals\';
format = '.tdms';
name = '20200819_171855_egg_probe2_01';%'20200721_141220_eggtest_DMA_01'%'20200529_162443_fish_egg_01'%20200604_174024_agarose_new_probe_1000kPa_close_01';

%%
%experimental data
pip_rad = 45e-6;%50e-6/2; %m
p_sens = 76.59;     %Pa/nm
Rc = 505e-6;        %m
%DMA
DMA_cycl =  [1 1 3 4 5];              %number of cycles
corr_time = 0;                    %lag before test starts
DMA_freq = [0.05 0.1 0.33 0.7 .85]; %testing frequencies
start_DMA = floor([22.43 43.95 55.04 66.43 74.35]*rate+corr_time); %start

%%
filename = strcat(path,name,format);        
dsens_struct = TDMS_getStruct(filename);

rate = 1000;
freq = 25;
Order = 3;
[b, a] = butter(Order, freq / (rate * 2), 'low');

ds_p = dsens_struct.Demodulated_data.pressure.data';
ds_pb = filter(b,a,ds_p);
ds_as = dsens_struct.Demodulated_data.aspirated_length.data';
ds_asb = filter(b,a,ds_as);
%%
demod_data = figure;
subplot(1,2,2)
plot(ds_asb/1.33,-ds_pb*p_sens)
xlabel('L_P [nm]')
ylabel('\DeltaP [Pa]')
%%
approx_mod = 3*1500*(100e-6/2)/((400e-9)*2*pi)*2.1/1e3; %kPa
legend('Butterworth')

time = 0.001:1/1000:length(ds_p)/1000;

subplot(1,2,1)
yyaxis right
plot(time,-ds_pb)
ylabel('\deltaOPL [nm]')
hold on
yyaxis left
plot(time,-ds_asb/1.33)
%plot(time_img,-1000*pos_img,'*')
hold off
xlabel('time [s]')
ylabel('\deltaOPL [nm]')
legend('displacement','pressure')
%legend('fiber-to-specimen','video tracking','pressure sensor')

%%
%find the ramp portion
%average diff(pressure) readout at rest
der = zeros(length(time),1);
der = filter(b,a,diff(-ds_pb*p_sens));
[minimum,min_der_idx] = min(der);
der = der(min_der_idx:end);
%min_der_idx = 0
pressure_start = mean(der(1:rate/5));
%start of aspiration
idx = find(der>-0.0001 & der<0.0001)%> 3.2*pressure_start);
idx_s = idx(1);
%find max(P)
pressure_max = max(der);
idx = find(der/pressure_max >= .9595);
idx_e = idx(1);
%plot the extrapolated aspiration curve
%  figure(2)
%  plot(ds_asb(idx_s+min_der_idx:idx_e+min_der_idx)/1.33,-p_sens*ds_pb(idx_s+min_der_idx:idx_e+min_der_idx))
%%
%calculate slope inf half space
corr = 2*pi/(3*2.1);
pf = polyfit(ds_asb(idx_s+min_der_idx:idx_e+min_der_idx)/1.33,...
    -p_sens*ds_pb(idx_s+min_der_idx:idx_e+min_der_idx),1);
EY_Theret = pf(1)*10^9/corr*pip_rad/1000 %kPa
%%
%slope with finite size correction
beta1 = 2.0142;
beta3 = 2.1187;
c_fit = beta1*(1-(pip_rad/Rc)^beta3)/3;
% pf2 = polyfit(ds_asb(idx_s+min_der_idx:idx_e+min_der_idx)/1.33,...
%     -p_sens*ds_pb(idx_s+min_der_idx:idx_e+min_der_idx),1);
pf2 = polyfit(ds_asb(2440:11830)/1.33,...
    -p_sens*ds_pb(2440:11830),1);
EY_Zhou = pf2(1)*10^9*pip_rad/c_fit/1000 %kPa
%plot over exp data and estimate the modulus
figure(demod_data)
subplot(1,2,2)
hold on
yfit = pf2(1)*ds_asb+pf2(2);
plot(ds_asb,yfit,'--r')
legend('Butterworth','linearized Zhou model');
%%
%test fit
%32270:42270
%44150:50760
%52660:60980
%62910:69870
%76020:84200
%60260:64500
%9250:29030
%29700:40120
%55780:60990
%start_DMA = [11840 106100 126600 136000 142400]+corr_time;
%start_DMA = [11080 29700 38860 50100 60190 66510]+corr_time;
%fin_DMA = [104600 122400 134100 141100 147000]+corr_time;
%fin_DMA = [27500 41070 46550 58380 65170 70500]+corr_time;

%DMA_freq = [0.05 0.1 0.35 0.751 1.01]; %testing frequencies
%DMA_freq = [0.01 0.1 0.35 0.75 1.0]; %testing frequencies

%DMA_cycl =  [1 2 2 3 4];
%DMA_cycl =  [1 2 3 4 5];
%start_DMA = floor([8.9 27.5 47.5 56.8 62.8]*rate+corr_time);   %timepoints
%start_DMA = floor([8.9 106 125 134 140.9]*rate+corr_time);   %timepoints
fin_DMA = ceil(start_DMA+DMA_cycl./DMA_freq*rate);
Es = zeros(length(DMA_freq),1);
El = zeros(length(DMA_freq),1);
for i=1:length(start_DMA)
    ds_asb_f = ds_asb(start_DMA(i):fin_DMA(i))/1.331;
    ds_pb_f = -ds_pb(start_DMA(i):fin_DMA(i))*p_sens;
    time_f = time(start_DMA(i):fin_DMA(i));                 
    [~,E1,E2,phi] = fitDynSweep(time_f',ds_asb_f/1e9,ds_pb_f,...
                    pip_rad,Rc,1,0.5,DMA_freq(i),'egg');
    Es(i) = E1;
    El(i) = E2;
end
%%
figure(5)
plot(DMA_freq,Es/1000,'o-k')
hold on
xlabel('frequency [Hz]')
plot(DMA_freq,abs(El)/1000,'*-k')
yyaxis left
ylabel('E''/E'''' [kPa]')
yyaxis right
plot(DMA_freq,abs(El)./Es,'.-r')
ylabel('tand\delta')
ax=gca;
ax.XScale = 'log';
ax.YColor = 'r'
grid on
legend E'' E'''' tan\delta
hold off
%%
% %Zhou model
% beta1 = 2.0142;
% beta2 = 2.1186; 
% beta3 = 2.1187;
% beta4 = -1.4409;
% beta5 = 0.3154;
% Rp = 150e-6; %m
% EY = 55000; %Pa
% ni = 0.5;
% G = EY/2*(1+ni);
% L = 0:0.01:2200;
% L = 1e-9*L; %m
% Rc = 305.6e-6; %m
% dP = G*(beta1*L/Rp+beta1*(L/Rp).^2).*(1-(Rp/Rc).^(beta3+beta4*L/Rp+beta5*(L/Rp).^2));
% figure(31)
% plot(L*1e9,dP);

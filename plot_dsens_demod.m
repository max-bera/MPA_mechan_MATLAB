close all
%%
%deltasens data
path = 'C:\Users\Massimiliano\Desktop\DeltaSens signals\';
format = '.tdms';
name = '20200819_171855_egg_probe2_01';
%%
%experimental data
pip_rad = 45e-6;    %m
p_sens = 76.59;     %Pa/nm
Rc = 505e-6;        %m
rate = 1000;        %Hz

%DMA
ncycles =  [1 1 3 4 5];    %number of cycles
%--------------these are also used for quasi-static test-------------------
dt = 2;                    %lag before test starts [s]
rampt = 10;                %time to reach the preload [s]
%--------------------------------------------------------------------------
creept = 10;               %creep time at target preload [s]
waitt = 2;                 %pause between oscillations [s]
freqs = [0.05 0.1 0.33 0.7 .85]; %testing frequencies
%%
%deltasens processing
filename = strcat(path,name,format);        
dsens_struct = TDMS_getStruct(filename);

freq = 25; %cutoff
Order = 3;
[b, a] = butter(Order, freq / (rate * 2), 'low');

ds_p = dsens_struct.Demodulated_data.pressure.data';
ds_pb = filter(b,a,ds_p);
ds_as = dsens_struct.Demodulated_data.aspirated_length.data';
ds_asb = filter(b,a,ds_as);
%%
%plot experimental data
demod_data = figure;
subplot(1,2,2)
plot(ds_asb/1.33,-ds_pb*p_sens)
xlabel('L_P [nm]')
ylabel('\DeltaP [Pa]')
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
%%
%calculate slope with finite size correction
LZs = ceil(dt*rate);
LZf = LZs+floor(rampt*rate);
EY_Zhou = fitLinZhou(ds_asb(LZs:LZf)/1.331,-p_sens*ds_pb(LZs:LZf),...
          pip_rad,Rc,1);
%%
%DMA fit
[Es,El,tand] = fitDMA(time',ds_asb/1.331/1e9,-ds_pb*p_sens,pip_rad,Rc,...
               freqs,ncycles,dt,rampt,creept,waitt,rate,'egg',1);

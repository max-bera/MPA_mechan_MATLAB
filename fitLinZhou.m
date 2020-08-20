function [EY_LZ,S] = fitLinZhou(Lp,dP,pip_rad,Rc,plotcheck)
%%
%fitLinZhou.m fits an aspiration experiment to a linearized version of the
%Zhou model
%
%INPUT (* are mandatory)
%Lp [m]        *    array of aspirated lengths
%dP [Pa]       *    array of corresponding pressure points
%pip_rad [m]   *    pipette radius
%Rc [m]             specimen radius (if empty assume infinite half space)
%plotcheck [0/1]    plot with experimetnal data and linear fit (optional)
%
%OUTPUT
%EY_LZ [Pa]         Elastic modulus
%S                  array of values for error estimation with polyval.m
%%
if nargin<3; help fitLinZhou; return; end
if nargin<4, Rc=100000000000000000000; end
if nargin<5, plotcheck = 0; end
%%
%slope with finite size correction
beta1 = 2.0142;
beta3 = 2.1187;
c_fit = beta1*(1-(pip_rad/Rc)^beta3)/3;
[pf2,S] = polyfit(Lp,dP,1);
EY_LZ = pf2(1)*pip_rad/c_fit; %Pa
%plot over exp data and estimate the modulus
if plotcheck == true
    fitLinZhou_vs_expdata = figure;
    figure(fitLinZhou_vs_expdata);
    plot(Lp,dP,'k')
    hold on
    yfit = pf2(1)*Lp+pf2(2);
    plot(Lp,yfit,'--r') 
    xlabel('L_p [m]');
    ylabel('\DeltaP [Pa]');
    legend('exp data','linearized Zhou model');
end
end
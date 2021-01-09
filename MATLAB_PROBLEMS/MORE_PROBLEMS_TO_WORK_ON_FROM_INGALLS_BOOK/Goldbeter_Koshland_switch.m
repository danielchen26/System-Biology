%file Goldbeter_Koshland_switch.m
%response of the zero-order ultrasensitive covalent modification cycle
%Figure 6.7

%determine curves for various K values

K1=1;
K2=1;
wstar=linspace(0,1,100)
vrat1=wstar.*(1-wstar+K1)./((1-wstar).*(wstar+K2))

K1=0.1;
K2=0.1;
wstar=linspace(0,1,100)
vrat2=wstar.*(1-wstar+K1)./((1-wstar).*(wstar+K2))

K1=0.01;
K2=0.01;
wstar001=linspace(0,1,1000)
vrat3=wstar001.*(1-wstar001+K1)./((1-wstar001).*(wstar001+K2))


%generate figure
figure(1)
set(gca,'fontsize',14)
plot(vrat1,wstar,'k',vrat2,wstar,'k',vrat3,wstar001, 'k','Linewidth', 3)
axis([0 3 0 1])
xlabel('Concentration of activating enzyme ([E_{1T}], arbitrary units)')
ylabel('Fractional abundance of active protein (w^*)')

%File exponential_decay.m
%Figure 2.4: Exponentially decaying concentration profiles

%time domain
t=linspace(0,10,100);

%concentration profiles
a1=exp(-t);
a2=2.*exp(-t);
a3=3.*exp(-t);

%generate plot
figure(1)
set(gca,'fontsize',14)
plot(t,a1, 'k', t,a2, 'k--', t,a3, 'k-.', 'LineWidth', 3)
axis([0,5,0,3.2])
xlabel('Time')
ylabel('Concentration')
legend('D=1', 'D=2', 'D=3')

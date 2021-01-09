%File gma_approximations.m
%Figure 3.13: Comparions of GMA and Michaelis-Menten rate laws

t=linspace(0,4);

figure(1)
set(gca,'fontsize',14)
hold on
plot(t, 2*t./(1+t), 'k', 'Linewidth', 3)
plot(t, t.^0.4, 'k--', 'Linewidth', 2)
xlabel('Substrate concentration (arbitrary units)')
ylabel('Reaction rate (arbitrary units)')
legend('Michaelis-Menten', 'GMA')
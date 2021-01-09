% file gprotein_bodeplots.m 
% Figure 6.21

%transfer function for system calculated separately (Maple code available
%upon request)

%transfer function:
num = [5 2858]*1.57142*10^18;
den = [125*10^8  125*5.7172*10^10+1*10^8  125*6.6687*10^9+1*5.7172*10^10  6.6687*10^9];

%define system:
sys=tf(num, den)

%produce plots
W=logspace(-4,1, 150);
[MAG,PHASE] = bode(sys,W);
MAG=squeeze(MAG);
PHASE=squeeze(PHASE);

figure(1)
set(gca,'fontsize',14)
loglog(W/(6.28318),MAG/1e9, 'k', 'Linewidth', 3);
%axis([1e-4 1e6 1e-6 10])
xlabel('Frequency (Hz)')
ylabel('Gain')
axis([1e-4/(6.28318) 1e1/(6.28318) 1e-2 1e3])

figure(2)
set(gca,'fontsize',14)
semilogx(W/(6.28318),PHASE, 'k', 'Linewidth', 3);
%axis([1e-4 1e6 1e-6 10])
xlabel('Frequency (Hz)')
ylabel('Phase Shift')
axis([1e-4/(6.28318) 1e1/(6.28318) -180 30])
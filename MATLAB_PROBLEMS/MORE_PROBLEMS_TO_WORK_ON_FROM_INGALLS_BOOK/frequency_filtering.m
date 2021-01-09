% file frequency_filtering.m
% low-pass and band-pass filters
% figure 6.23

function frequncy_filtering

clear all

%declare variables
global y
global N
global Tend
global m;

%generate input signal
N=1000;
Tend=30;

y=zeros(N,1);
t=linspace(0,Tend,N);

for i=1:N
    y(i)=sin(t(i))+ 1*randn;
end



figure(1)
set(gcf,'Position',[0 0 900 800])

%generate Figure 6.22A
subplot(3,3,2)
set(gca,'fontsize',12)
plot(t,y, 'k', 'Linewidth', 2)
xlabel('Time (s)')
axis([0 30 -5 5])
h=text(-10,-1.5, 'input', 'Fontsize', 18);
set(h, 'rotation', 90)


%simulate low-pass filter in Figure 6.22B
x0=0;
m=-10;
ODEFUN=@lowfilterddt;
[tt,ss]=ode45(ODEFUN, [0,Tend], x0);

%generate figure 6.22C
subplot(3,3,7)
set(gca,'fontsize',12)
plot(tt,10*ss, 'k', 'Linewidth', 2)
xlabel('Time (s)')
axis([0 30 -2 2])
h=text(-10,-0.75, 'output', 'Fontsize', 18);
set(h, 'rotation', 90)

%simulate low-pass filter in Figure 6.22D
x0=0;
m=-0.1;
ODEFUN=@lowfilterddt;
[tt,ss]=ode45(ODEFUN, [0,Tend], x0);

%generate figure 6.22E
subplot(3,3,8)
set(gca,'fontsize',12)
plot(tt,ss, 'k', 'Linewidth', 3)
axis([0 30 -2 2])
xlabel('Time (s)')

%simulate band-pass filter in Figure 6.22F
x0=[0 0]';
ODEFUN=@bandfilterddt;
[tt,ss]=ode45(ODEFUN, [0,Tend], x0);

%generate figure 6.22G
subplot(3,3,9)
set(gca,'fontsize',12)
plot(tt,ss(:,2), 'k', 'Linewidth', 2)
axis([0 30 -4 4])
xlabel('Time (s)')



%Bode plots

%generate magnitude Bode plot for low-pass filter
W=logspace(-4,6, 150);
num =10;
m = 10;
den=[1 m];
sys=tf(num,den);
[MAG,PHASE] = bode(sys,W);
MAG=squeeze(MAG);

%generate Figure 6.22B
subplot(3,3,4)
set(gca,'fontsize',12)
loglog(W/(6.28318),MAG, 'k', 'Linewidth', 3);
axis([1e-4/(6.28318) 1e6/(6.28318) 1e-6 10])
xlabel('Frequency (Hz)')
ylabel('Gain')
set(gca, 'Xtick', [1e-3 1e0 1e3])
h2=text(1e-10,1e-4, 'system', 'Fontsize', 18);
set(h2, 'rotation', 90)


%generate magnitude Bode plot for low-pass filter
num =0.1;
m = 0.1;
den=[1 m];
sys=tf(num,den);
[MAG,PHASE] = bode(sys,W);
MAG=squeeze(MAG);

%generate Figure 6.22D
subplot(3,3,5)
set(gca,'fontsize',12)
loglog(W/(6.28318),MAG, 'k', 'Linewidth', 3);
axis([1e-4/(6.28318) 1e6/(6.28318) 1e-6 10])
xlabel('Frequency (Hz)')
ylabel('Gain')
set(gca, 'Xtick', [1e-3 1e0 1e3])

%generate magnitude Bode plot for band-pass filter
A=[-100/1 0; 100/1 -10000/5];
B=[1000/1 -1000/1]';
C=[0 1];
D=0;
[num den]=ss2tf(A,B,C,D);
sys=tf(num,den);
[MAG,PHASE] = bode(sys,W);
MAG=squeeze(MAG);

%generate Figure 6.22F
subplot(3,3,6)
set(gca,'fontsize',12)
loglog(W/(6.28318),MAG, 'k', 'Linewidth', 3);
axis([1e-4/(6.28318) 1e6/(6.28318) 1e-6 10])
xlabel('Frequency (Hz)')
ylabel('Gain')
set(gca, 'Xtick', [1e-3 1e0 1e3])

end

%dynamics for low-pass filtering
function dS=lowfilterddt(t,x)

global y
global N
global Tend
global m

dS = m*x + y(max(ceil((t/Tend)*N),1));

end



%dynamics for band-pass filtering
function dS=bandfilterddt(t,x)

global y
global N
global Tend

dS =[-100*x(1) + 1000*y(max(ceil((t/Tend)*N),1)), 100*x(1) - 10000/(5)*x(2) - 1000*y(max(ceil((t/Tend)*N),1))]';

end



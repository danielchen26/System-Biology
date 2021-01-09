% file resonant_system.m
% filtering behaviour of a resonant system
% figure 6.22


function resonant_system

clear all

%delcare variables
global y
global N
global Tend


%Linear system with matrices as follows:
% A =
% 
%      0   -10
%      1     0
% 
% 
% B =
% 
%      1
%      0
% 
% 
% C =
% 
%      0     1
% 
% 
% D =
% 
%      0



%generate white noise signal
N=3000;
Tend=340;

y=zeros(N,1);
t=linspace(0,Tend,N);

randn('state',1)
for i=1:N
    y(i)=randn;
end
for i=1:N
    y(i)=randn;
end

%run simulation
x0=[0 0]';
ODEFUN=@resonantddt;
[tt,ss]=ode45(ODEFUN, [0,Tend], x0);

%generare figure 6.22A
figure(1)
set(gca,'fontsize',18)
plot(t-300,y, 'k', 'Linewidth', 2)
axis([0 30 -4 4])
xlabel('Time (s)')
ylabel('Input')
str1(1) = {'A'};
text(-3.5,4,str1, 'Fontsize', 40)

%generate figure 6.22C
figure(3)
set(gca,'fontsize',18)
plot(tt-300,ss(:,2), 'k', 'Linewidth', 3)
axis([0 30 -3 3])
xlabel('Time (s)')
ylabel('Output')
str1(1) = {'C'};
text(-3.5,3,str1, 'Fontsize', 40)



%generate Bode magnitude plot
W=logspace(-1,2, 1500);
num=1;
den=[1 0 10];
sys=tf(num,den);
[MAG,PHASE] = bode(sys,W);
MAG=squeeze(MAG);

%generate figure 6.22B
figure(2)
set(gca,'fontsize',18)
loglog(W/(6.28318),MAG, 'k', 'Linewidth',3);
axis([1e-1/(6.28318) 1e2/(6.28318) 1e-3 1e2])
xlabel('Frequency (Hz)')
ylabel('Gain')
str1(1) = {'B'};
text(0.006,100,str1, 'Fontsize', 40)

end

%dynamics
function dS=resonantddt(t,x)

global y
global N
global Tend

dS =[-10*x(2) + y(max(ceil((t/Tend)*N),1)), x(1)]';

end


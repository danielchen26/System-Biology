%file synchronized_oscillators.m
%Hasty synthetic oscillator model
%from Hasty et al. (2002) Phys. Rev. Lett. 88 148101
%Figure 7.31

function synchronized_oscillators

clear all

%declare model parameters
global alpha
global sigma
global gammax
global gammay
global ay
global D

%assign parameter values
alpha=11
sigma=2
gammax=0.2
gammay=0.012
ay=0.2
D=0.015

%set simulation parameters
ODEFUN=@syncddt;
options=odeset('Refine', 6);
Tend=500;

%set initial condition. State is [x1 y1 x2 y2] 
x0=[0.3963    2.3346  0.5578    1.9317]';

%generate simulation
[t,s]=ode45(ODEFUN, [0 Tend], x0, options);

%produce figure 7.31
  figure(1)
 set(gca,'fontsize',14)
c=[0.6 0.6 0.6; 0.6 0.6 0.6; 0 0 0;0 0 0]
set(0,'DefaultAxesColorOrder',c)
 plot(t,s(:,3), t,s(:,4), '--', t,s(:,1), t,s(:,2),'--','Linewidth', 3)
 xlabel('Time (arbitrary units)')
 ylabel('Concentration (arbitrary units)')
 legend('X_1', 'Y_1', 'X_2', 'Y_2')
 

end

%dynamics
function dS=syncddt(t,s)

global alpha
global sigma
global gammax
global gammay
global ay
global D

x=s(1);
y=s(2);
xx=s(3);
yy=s(4);

%original model
xddt= (1+x^2+alpha*sigma*x^4)/((1+x^2+sigma*x^4)*(1+y^4)) - gammax*x + D*(xx-x);
yddt= (ay)*((1+x^2+alpha*sigma*x^4)/((1+x^2+sigma*x^4)*(1+y^4))) - gammay*y;

xxddt= (1+xx^2+alpha*sigma*xx^4)/((1+xx^2+sigma*xx^4)*(1+yy^4)) - gammax*xx + D*(x-xx);
yyddt= (ay)*((1+xx^2+alpha*sigma*xx^4)/((1+xx^2+sigma*xx^4)*(1+yy^4))) - gammay*yy;

dS =[xddt,  yddt, xxddt, yyddt]';

end


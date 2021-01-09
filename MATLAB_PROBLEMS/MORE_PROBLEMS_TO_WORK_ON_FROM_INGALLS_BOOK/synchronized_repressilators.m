% file synchronized_repressilators.m
% model of synchronized repressilator networks
% from Garcia-Ojalvo et al. (2004) PNAS 101 pp. 10955-10960
% problem 7.8.19


function synchronized_repressilators

clear all

%declare model parameters
global alpha
global alpha0
global n
global beta
global kappa

global ks0;
global ks1;
global eta
global Q

%assign parameter values
alpha0=0; 
alpha=216;  
n=2;
beta=1;

kappa=20;
ks0=1;
eta=2;
ks1=0.01;
Q=0.9;

%set simulation parameters
ODEFUN=@repressilatorddt;
options=odeset('Refine', 6);
Tend=100;

%set initial condition, state is [m11 p11 m21 p21 m31 p31 s1 m12 p12 m22 p22
%m32 p32 s2]
x0=[0  10 0 0 0 0 0 0  0 0 10 0 0 0]';

%run simluation and plot
[t,s]=ode45(ODEFUN, [0,Tend], x0, options);
figure(1)
plot(t,s(:,2), t,s(:,9),'Linewidth', 3)

end

%dynamics
function dS=repressilatorddt(t,s)


global alpha
global alpha0
global n
global beta
global kappa

global ks0;
global ks1;
global eta
global Q



ma1=s(1);
pa1=s(2);
mb1=s(3);
pb1=s(4);
mc1=s(5);
pc1=s(6);
s1=s(7);
ma2=s(8);
pa2=s(9);
mb2=s(10);
pb2=s(11);
mc2=s(12);
pc2=s(13);
s2=s(14);


se=Q*((1/2)*(s1+s2));

ma1dt= alpha0 + alpha/(1+pc1^n) - ma1;
pa1dt=beta*(ma1-pa1);
mb1dt= alpha0 + alpha/(1+pa1^n) - mb1;
pb1dt=beta*(mb1-pb1);
mc1dt= alpha0 + alpha/(1+pb1^n) + kappa*s1/(1+s1) - mc1;
pc1dt=beta*(mc1-pc1);
s1dt=-ks0*s1+ks1*pa1-eta*(s1-se);

ma2dt= alpha0 + alpha/(1+pc2^n) - ma2;
pa2dt=beta*(ma2-pa2);
mb2dt= alpha0 + alpha/(1+pa2^n) - mb2;
pb2dt=beta*(mb2-pb2);
mc2dt= alpha0 + alpha/(1+pb2^n) + kappa*s2/(1+s2) - mc2;
pc2dt=beta*(mc2-pc2);
s2dt=-ks0*s2+ks1*pa2-eta*(s2-se);

dS =[ma1dt,  pa1dt, mb1dt, pb1dt, mc1dt, pc1dt, s1dt, ma2dt,  pa2dt, mb2dt, pb2dt, mc2dt, pc2dt, s2dt]';

end


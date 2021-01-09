% file population_control.m
% model of population control network
% from You et al. (2004) Nature 428 pp. 868-871 
% problem 7.8.20 (a)


function population_control

%declare model parameters
global k;
global Nm;
global d;
global ke;
global de;
global va;
global da;
 
%assign parameter values
k=0.97;
Nm=1.24e9;
d=0.004;
ke=5;
de=2;
va=4.8e-7;
da=0.64;

%set simulation parameters
ODEFUN=@popddt;
Tend=100;

%set initial condition. state: [N, E, A]
S0=[1 1 1];


[t,S]=ode45(ODEFUN, [0,Tend], S0);


figure(1)
set(gca, 'fontsize', 14)
plot(t, S(:,1))
xlabel('Time (h)')
    

end

%dynamics
function dS = popddt(t,S)

global k;
global Nm;
global d;
global ke;
global de;
global va;
global da;
 
N=S(1);
E=S(2);
A=S(3);


dNdt=k*N*(1-N/Nm)-d*E*N;
dEdt=ke*A-de*E;
dAdt=va*N-da*A;

dS=[dNdt;dEdt;dAdt];

end
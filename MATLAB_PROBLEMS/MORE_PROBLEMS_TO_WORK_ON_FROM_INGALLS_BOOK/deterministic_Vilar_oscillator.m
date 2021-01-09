% file deterministic_vilar_oscillator.m
% deterministic simulation of genetic oscillator
% from Vilar et al. (2002) PNAS 99 5988-5992.
% problem 7.8.27

function vilar
clear all

%declare parameter values
global ga
global ba
global Ka
global alpha0
global da
global kc
global gr
global br
global Kr
global dr

%assign parameter values
eps=0.1;
ga=250; 
ba=5;
Ka=0.5; 
alpha0=0.1;
da=1; 
kc=200;
gr=50; 
br=10;
Kr=1; 
dr=eps*da;

%set simulation parameters
ODEFUN=@vilarddt;
options=odeset('Refine', 6);
Tend=50;

%set initial condition: state [A R C]
x0=[5  10  35]';

%run simulation
[t,s]=ode23s(ODEFUN, [0 Tend], x0, options);

%plot figure
 figure(3)
 plot(t,s(:,1),'k', t,s(:,2), 'k--',t,s(:,3),'g', 'Linewidth', 3)
 xlabel('Time (min)')
 ylabel('Concentration (arbitrary units)')
 legend('A', 'R', 'C')



end

%dynamics
function dS=vilarddt(t,s)

global ga
global ba
global Ka
global alpha0
global da
global kc
global gr
global br
global Kr
global dr


A=s(1);
R=s(2);
C=s(3);


Addt=ba*(ga/ba)*(alpha0+(A/Ka))/(1+(A/Ka)) - da*A - kc*A*R;
Rddt=br*(gr/br)*((A/Kr))/(1+(A/Kr)) - dr*R - kc*A*R + da*C;
Cddt= kc*A*R - da*C;
            

dS =[Addt,  Rddt, Cddt]';

end


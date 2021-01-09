% file predator_prey.m
% model of predator-prey network
% from Balagadde et al. (2008) Molecular Systems Biology 4
% problem 7.8.20 (b)


function predatorprey

%declare model parameters
global k1;
global k2;
global cm;
global d1;
global d2;
global K1;
global beta;
global D;
global gamma1;
global gamma2;
global del1;
global del2;

%assign parameter values
 k2=0.4;
 cm=1e5;
 d1=1; 
 K1=10;
 beta=2;
 D=0.2; 
 gamma1=0.1;
 gamma2=0.1; 
 del1=0.017;
 del2=0.11;
 k1=0.8;
 d2=0.3;
 
 
Tend=500;

%state: [c1 c2 a1 a2]
S0=[10000 10000 1 1];

ODEFUN=@predpreyddt;
[t,S]=ode45(ODEFUN, [0,Tend], S0);

figure(1)
set(gca, 'fontsize', 14)
plot(t, S(:,1), 'k', t, S(:,2),'k-.', 'LineWidth', 3)
legend('c1', 'c2')
xlabel('Time (h)')
   
end

%dynamics
function dS = predpreyddt(t,S)

global k1;
global k2;
global cm;
global d1;
global d2;
global K1;
global beta;
global D;
global gamma1;
global gamma2;
global del1;
global del2;

c1=S(1);
c2=S(2);
a1=S(3);
a2=S(4);

dc1dt=k1*c1*(1-(c1+c2)/cm)-d1*c1*K1/(K1+a2^beta) - D*c1;
dc2dt=k2*c2*(1-(c1+c2)/cm)-d2*c2*a1^beta/(K1+a1^beta) - D*c2;
da1dt=gamma1*c1-(del1+D)*a1;
da2dt=gamma2*c2-(del2+D)*a2;

dS=[dc1dt;dc2dt;da1dt;da2dt];

end
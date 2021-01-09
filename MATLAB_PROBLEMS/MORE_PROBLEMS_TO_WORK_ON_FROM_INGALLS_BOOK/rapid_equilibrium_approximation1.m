%File rapid_equilibrium_approximation1.m
%Figure 2.11 (Fig 1) & Figure 2.12 (Fig 2) 
%Simulation and rapid equilibrium approximation

function rapid_equilibrium_approximation1
%declare model parameters
global k1;
global km1;
global k2;
%assign parameter values
k1=9;
km1=12;
k2=2;
%set final time for simulation
Tend=3;
%set initial conditions for original model: S(1)=a, S(2)=b
S0=[0,10];
%define right-hand side (function defined below)
ODEFUN=@dSdt_original;
%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], S0);
%plot figure
figure(1)
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', t, S(:,2), 'k--', 'LineWidth',3)
xlabel('Time (arbitrary units)')
ylabel('Concentration (arbitrary units)')
legend('a', 'b')

%set initial condition for approximate model: S(1)=c
S0=[10];
%define right-hand side (function defined below)
ODEFUN=@dSdt_approximation;
%run simulation
[t_approx,S_approx]=ode45(ODEFUN, [0,Tend], S0);
%plot figure (a_tilde and b_tilde are determined in the plot function)
figure(2)
set(gca,'fontsize',14)
hold on
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', 'LineWidth',2)
plot(t, S(:,2), 'k--', 'LineWidth',2)
plot(t_approx, km1/(k1+km1)*S_approx(:,1), 'k-.', t_approx, k1/(k1+km1)*S_approx(:,1),'k:', 'LineWidth',2)
legend('a (original model)', 'b (original model)', 'a (reduced model)', 'b (reduced model)')
xlabel('Time (arbitrary units)')
ylabel('Concentration (arbitrary units)')

end

%declare right-hand-side for original model
function dS = dSdt_original(t,S)

global k1;
global km1;
global k2;

a=S(1);
b=S(2);

dS=[-k1*a+km1*b, k1*a-km1*b-k2*b]';

end

%declare right-hand-side for approximate model
function dS = dSdt_approximation(t,S)

global k1;
global km1;
global k2;

c=S(1);

dS=[-k2*k1/(k1+km1)*c]';

end
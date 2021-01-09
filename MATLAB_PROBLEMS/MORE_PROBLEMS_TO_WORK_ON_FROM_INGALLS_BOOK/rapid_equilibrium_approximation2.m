%File rapid_equilibrium_approximation2.m
%Figure 2.13 Rapid equilibrium approximation

function rapid_equilibrium_approximation2
%declare model parameters
global k0;
global k1;
global km1;
global k2;
%assign parameter values
k0=5;
k1=20;
km1=12;
k2=2;
%set final time for simulation
Tend=3;
%set initial conditions for original model: S(1)=a, S(2)=b
S0=[8,4];
%define right-hand side (function defined below)
ODEFUN=@dSdt_original;
%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], S0);

%set initial condition for approximate model: c
S0=[12];
%define right-hand side (function defined below)
ODEFUN=@dSdt_approximate;
%run simulation
[t_approx,S_approx]=ode45(ODEFUN, [0,Tend], S0);
%calcuate approximate concentrations
a_tilde=km1/(km1+k1).*S_approx;
b_tilde=k1/(km1+k1).*S_approx;
%plot figure
figure(1)
hold on
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', 'LineWidth',2)
plot(t, S(:,2), 'k--', 'LineWidth',2)
plot(t_approx, a_tilde, 'k-.', 'LineWidth',2)
plot(t_approx, b_tilde, 'k:', 'LineWidth',2)
legend('a (original model)', 'b (original model)', 'a (reduced model)', 'b (reduced model)')
xlabel('Time (arbitrary units)')
ylabel('Concentration (arbitrary units)')


end


%declare right-hand-side for original model
function dS = dSdt_original(t,S)

global k0;
global k1;
global km1;
global k2;

%a=S(1)
%b=S(2)

dS=[k0 - k1*S(1) + km1*S(2), k1*S(1) - km1*S(2) - k2*S(2)]';

end

%declare right-hand-side for approximate model
function dR = dSdt_approximate(t,S)

global k0;
global k1;
global km1;
global k2;

%c=S(1)

dR=[k0 - (k2*k1)/(km1+k1)*S(1)];

end

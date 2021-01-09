%File quasi_steady_state_approximation.m
%Figure 2.14 Quasi-steady-state approximation

function quasi_steady_state_approximation
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
Tend=4;
%set initial condition for original model [a,b]
S0=[8,4];
%set right-hand-side for original model
ODEFUN=@dSdt_original;
%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], S0);

%set initial condition for approximate model [b_tilde]
S0=[235/32];
%set right-hand-side for approximate model
ODEFUN=@dSdt_approximate;
%run simulation
[t_approx,S_approx]=ode45(ODEFUN, [0,Tend], S0);
%determine approximate concentration a_tilde
a_tilde=(k0+km1.*S_approx)/k1;

%generate plot
figure(1)
hold on
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', 'Linewidth', 3)
plot(t, S(:,2), 'k--', 'Linewidth', 3)
plot(t_approx, a_tilde(:,1), 'k-.', t_approx, S_approx, 'k:', 'Linewidth', 3)
legend('a (original model)', 'b (original model)', 'a (reduced model)', 'b (reduced model)')
xlabel('Time (arbitrary units)')
ylabel('Concentration (arbitrary units)')
end

%dynamics of full model
function dS = dSdt_original(t,S)

global k0;
global k1;
global km1;
global k2;

%S(1)=a
%S(2)=b

dS=[k0 - k1*S(1) + km1*S(2), k1*S(1) - km1*S(2) - k2*S(2)]';

end

%reduced dynamics for qssa
function dS = dSdt_approximate(t,S)

global k0;
global k1;
global km1;
global k2;

%S(1)=b_tilde

dS=[k0 - k2*S(1)];

end






function Problem2_4_6
%declare model parameters
global k;
%assign parameter values
k=1;
%set final time for simulation
Tend=10;
%set initial condition 
S0=[0];
%set right-hand-side for original model
ODEFUN=@dSdt;
%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], S0);

%generate plot
figure(1)
hold on
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', 'Linewidth', 3)
xlabel('Time (arbitrary units)')
end

function dS = dSdt(t,S)

global k

dS=[k*(-S(1)+1)]';

end
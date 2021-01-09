%File network_example.m
%Figure 2.9 Numerical Simulation of network

function network_example

%set initial condition
S0=[0,0,0,0];
%declare right-hand-side function (defined below)
ODEFUN=@netexddt;
%perform simulation over time-interval [0,10]
[t,S]=ode45(ODEFUN, [0,10], S0);

%generate plot
figure(1);
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', t, S(:,2), 'k--',t, S(:,3), 'k:',t, S(:,4), 'k-.', 'Linewidth', 3)
legend('A', 'B', 'C', 'D')
axis([0 4 0 1])
xlabel('Time (sec)')
ylabel('Concentration (mM)')
 
end

%define right-hand-side of differential equations
function dS=netexddt(t,S)
%a=S(1)
%b=S(2)
%c=S(3)
%d=S(4)
dS=[3 - 2*S(1) - 2.5*S(1)*S(2),...
	   2*S(1) - 2.5*S(1)*S(2),...
       2.5*S(1)*S(2) - 3*S(3),...
       2.5*S(1)*S(2) - 4*S(4)]';
end
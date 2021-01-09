%File problem3_7_5.m
%Problem 3.7.5 part (a)

function chain

%declare enzyme total as a model parameter
global v0;
global Vm1;
global Vm2;
global Vm3;
global Km1;
global Km2;
global Km3;


%assign values to model parameters
v0=2;
Vm1=9;
Vm2=12;
Vm3=15;
Km1=1;
Km2=0.4;
Km3=3;


%set initial condition. State vector is [s1, s2, s3]
S0=[0.3,0.2,0.1];

%set right-hand-side
ODEFUN=@chainddt;
%run simulation
[t,S]=ode45(ODEFUN, [0,2], S0);

%generate plot
figure(1);
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', t, S(:,2), 'k--',t, S(:,3), 'k:', 'Linewidth', 3)
legend('s1', 's2', 's3')
     axis([0 2 0 0.8])
     xlabel('Time (arbitrary units)')
     ylabel('Concentration (arbitrary units)')
     

     

end

%dynamics of reaction chain
function dS=chainddt(t,S)

global v0;
global Vm1;
global Vm2;
global Vm3;
global Km1;
global Km2;
global Km3;


     dS=[v0 - Vm1*S(1)/(Km1+S(1)), 
	   Vm1*S(1)/(Km1+S(1))-Vm2*S(2)/(Km2+S(2)), 
	   Vm2*S(2)/(Km2+S(2))-Vm3*S(3)/(Km3+S(3))];
end



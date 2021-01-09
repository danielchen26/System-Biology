%File michaelis_menten.m
%Figure 3.3 Michaelis-Menten kinetics

function MichaelisMenten

%declare enzyme total as a model parameter
global et;
global k1;
global km1;
global k2;

%assign values to model parameters
k1=30;
km1=1;
k2=10;
et=1;


%%%Original model simulation

%set initial condition. State vector is [s, c, p]
S0=[5,0,0];

%set right-hand-side
ODEFUN=@MMdynddt;
%run simulation
[t,S]=ode45(ODEFUN, [0,1], S0);

%generate plot
figure(1);
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', t, S(:,2), 'k--',t, S(:,3), 'k:',t, et-S(:,2), 'k-.', 'Linewidth', 3)
legend('s', 'c', 'p', 'e')
     axis([0 1 0 5])
     xlabel('Time (arbitrary units)')
     ylabel('Concentration (arbitrary units)')
     
str1(1) = {'A'};
text(-.12,5,str1, 'Fontsize', 40)
     

     
%%%%%%Reduced model
     
%set initial condition. State is s
Sred0=[5]

%set right-hand-side
ODEFUN=@redMMdynddt;
%run simulation
[t2,S2]=ode45(ODEFUN, [0,2], Sred0);

%generate plot
figure(2)
set(gca,'fontsize',14)
plot(t, S(:,1), 'k', t, S(:,3), 'k:',t2, S2(:), 'k--',t2, 5-S2(:), 'k-.', 'Linewidth', 3)
legend('s (full model)', 'p (full model)', 's (reduced model)', 'p (reduced model)')
     axis([0 1 0 5])
     xlabel('Time (arbitrary units)')
     ylabel('Concentration (arbitrary units)')
     
     str1(1) = {'B'};
text(-.12,5,str1, 'Fontsize', 40)
 

end

%dynamics of full model
function dS=MMdynddt(t,S)

global et;
global k1;
global km1;
global k2;

     dS=[-k1*S(1)*(et-S(2)) + km1*S(2), 
	   -km1*S(2) + k1*S(1)*(et-S(2))-k2*S(2), 
	   k2*S(2)];
end


%dynamics of reduced model
function dS=redMMdynddt(t,S)

global et;
global k1;
global km1;
global k2;


dS=[-k1*k2*et*S/(km1+k2+k1*S)];
end

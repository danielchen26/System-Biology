% file NfkB_signaling.m
% model of Nf-kB signaling pathway
% from Krishna et al. (2006) PNAS 103, pp. 10840-10845
% problem 7.8.15

function NFkB_signaling

%declare model parameters
global A;
global eps;
global B;
global del;
global C;

%assign parameter values
A=0.007;
B=954.5;
C=0.035;
del=0.029;
eps=2e-5;

%set simulation parameters
ODEFUN=@NFkBddt;
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Refine', 3);
Tend=20;

%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], [0.0014    0.0207    0.0127]);

%produce plot
figure(1)
plot(t, S(:,1), t, S(:,2),t, S(:,3), 'k', 'LineWidth',3)
xlabel('Time')
ylabel('Concentration')
legend('Nn', 'Im', 'I')


end




%dynamics
function dS = NFkBddt(t,S)

global A;
global eps;
global B;
global del;
global C;

Nn=S(1);
Im=S(2);
I=S(3);

dNndt=A*(1-Nn)/(eps+I) - B*I*Nn/(del+Nn);
dImdt=Nn^2-Im;
dIdt=Im-C*(1-Nn)*I/(eps+I);

dS=[dNndt dImdt dIdt]';

    
end



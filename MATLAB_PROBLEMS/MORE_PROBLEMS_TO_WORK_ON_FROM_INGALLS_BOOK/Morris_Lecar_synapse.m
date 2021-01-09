%file morris_lecar_synapse.m
%Morris-Lecar model of synapse onto excitable barnacle muscle fiber
%adapted from Morris and Lecar (1981) Biophysical Journal 35 pp. 193-231
%Figure 8.9

function morris_lecar_synapse

%declare model parameters
global C;
global gbarCa;
global ECa;
global gbarK;
global EK;
global gleak;
global Eleak;
global v1;
global v2;
global v3;
global v4;
global phi;
global Iapplied;
global gsyn;
global Esyn;
global alpha;
global beta;
global tau;
global Tmax;
global Vsynstar;


%assign parameter values
C= 20 ; %microfarad/cm^2 
gbarCa=4.4; % millisiemens/ cm^2 
ECa=120; %millivolts
gbarK=8;% millisiemens/ cm^2 
EK=-84; %millivolts
gleak=2;% millisiemens/ cm^2 
Eleak=-60;%millivolts
v1=-1.2; %millivolts
v2= 18 ; %millivolts
v3= 2 ; %millivolts
v4= 30; %millivolts
phi = 0.04 % per millisecond
tau=0.8;

Iapplied=0;

gsyn=1; % millisiemens/ cm^2 
Esyn= 150; % millisiemens/ cm^2 
alpha=1 %/ms
beta=0.12 %/ms
Tmax=1;
Vsynstar=27% mV

%simulation time [0, Tend]
Tend=200;

%set initial conditions.  State is (V,w)
%rest:
S00=[-60.8554    0.0149  0  -60.8554    0.0149];
%super-threshold
S10=[ -10  0.0149 0  -60.8554    0.0149];

%run simulation
ODEFUN=@synapse_ddt;
[t1,S1]=ode15s(ODEFUN, [-10 100], S00);
[t2,S2]=ode15s(ODEFUN, [t1(length(t1)) 200 + t1(length(t1))], S10);

%produce Figure 8.9A
figure(10)
subplot(2,1,1)
set(gca, 'fontsize', 14)
plot([t1;t2], [S1(:,1); S2(:,1)],'k',[t1;t2], [S1(:,4); S2(:,4)],'k--','Linewidth', 3)
xlabel('Time (msec)', 'fontsize',12)
ylabel('Membrane voltage (mV)', 'fontsize',12)
legend('V_{pre}', 'V_{post}');
axis([0 300 -75 40])
str1(1) = {'A'};
text(-50,40,str1, 'Fontsize', 40)

%produce Figure 8.9B
subplot(2,1,2)
set(gca, 'fontsize', 14)
plot([t1;t2], [S1(:,3); S2(:,3)],'k', 'Linewidth', 3)
xlabel('Time (msec)', 'fontsize',12)
ylabel('Current (pA/cm^2)', 'fontsize',12)
legend('I_{syn}');
axis([0 300 -0.1 1])
str1(1) = {'B'};
text(-50,0.9,str1, 'Fontsize', 40)

end



%dynamics
function dS = synapse_ddt(t,S)

global C;
global gbarCa;
global ECa;
global gbarK;
global EK;
global gleak;
global Eleak;
global v1;
global v2;
global v3;
global v4;
global phi;
global Iapplied;
global gsyn;
global Esyn;
global alpha;
global beta;
global tau;
global Tmax;
global Vsynstar;


%locally define state variables:
Vpre=S(1);
wpre=S(2);
ssyn=S(3);
Vpost=S(4);
wpost=S(5);

%local functions:
mpre_inf = (0.5*(1+tanh((Vpre-v1)/v2)));
wpre_inf = (0.5*(1+tanh((Vpre-v3)/v4)));

mpost_inf = (0.5*(1+tanh((Vpost-v1)/v2)));
wpost_inf = (0.5*(1+tanh((Vpost-v3)/v4)));

if Vpre>=Vsynstar
    T=Tmax;
else
    T=0;
end

ddt_Vpre = (1/C)*(gbarCa*mpre_inf*(ECa-Vpre) + gbarK*wpre*(EK-Vpre) + gleak*(Eleak-Vpre));
ddt_wpre = phi*(wpre_inf-wpre)/(tau);
ddt_ssyn = alpha*T*(1-ssyn) - beta*ssyn;
ddt_Vpost = (1/C)*(gbarCa*mpost_inf*(ECa-Vpost) + gbarK*wpost*(EK-Vpost) + gleak*(Eleak-Vpost)+gsyn*ssyn*(Esyn-Vpost));
ddt_wpost = phi*(wpost_inf-wpost)/(tau);


dS=[ddt_Vpre; ddt_wpre; ddt_ssyn; ddt_Vpost; ddt_wpost];

end
 

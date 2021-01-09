%file arabidopsis_branch.m
%from Curien et al., 2003, Eur. J. Biochemistry 270, p. 4615-4627
%Problem 5.6.7

function arabidopsis_branch

%declare parameters
global kcatCGS;
global KmCGScys;
global KmCGSPhser;
global KiCGSPi;

global KiTSPi;
global JPhser;

global Cys;
global CGS;
global Pi;
global ADOMET;
global TS;

%assign parameter values

kcatCGS=30; %/sec
KmCGScys=460; %muM
KmCGSPhser=2500; %muM
KiCGSPi=2000; %muM
 
KiTSPi =1000; %muM
JPhser=0.3; %muM/sec
 
Cys=250;%muM
CGS=0.7;%muM
Pi=10000;%muM
ADOMET=20; %muM
TS=5; %muM


%set up simulation options
ODEFUN=@branchddt;
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Refine', 3);
%set simulation time
Tend=5000;

%run simulation
[t,S]=ode45(ODEFUN, [0,Tend], [150]);

%plot results
figure(21)
plot(t, S, 'k', 'LineWidth',3)
%axis([0 Tend 0 5])
xlabel('Time')
ylabel('Concentration')

end


%dynamics for bistable metabolic chain
function dS = branchddt(t,S)

global kcatCGS;
global KmCGScys;
global KmCGSPhser;
global KiCGSPi;

global KiTSPi;
global JPhser;
global Cys;
global CGS;
global Pi;
global ADOMET;
global TS;

Phser=S;

kcatappCGS=kcatCGS/(1+KmCGScys/Cys);

KmCGSapp=(KmCGSPhser/(1+KmCGScys/Cys))*(1+Pi/KiCGSPi);

kTS=5.4*1e-5+ 6.2*1e-3*ADOMET^2.9/(32^2.9+ADOMET^2.9);

 vThr=TS*kTS/(1+Pi/KiTSPi)*Phser;
 
 vcystathionine=kcatappCGS*CGS*Phser/(KmCGSapp+Phser);

dS=[JPhser - vThr-vcystathionine];

    
end



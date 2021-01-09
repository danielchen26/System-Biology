% file mapk_pathway.m
% model of MAPK signalling pathway
% from Kholodenkho (2000) Euro. J. Biochem. 267 pp 1583-1588
% problem 6.8.8

function mapk_pathway
%from Kholodenkho, 2000

%declare model parameters
global V1;
global KA;
global KI;
global KM1;
global V2;
global KM2;
global kcat3;
global KM3;
global kcat4;
global KM4;
global V5;
global KM5;
global V6;
global KM6;
global kcat7;
global KM7;
global kcat8;
global KM8;
global V9;
global KM9;
global V10;
global KM10;


%assign model parameters
 
 KA=0;
 KI=0;
 
 V1=2.5;
 KM1=10;
 V2=0.25;
 KM2=8;
 kcat3=0.025;
 KM3=15;
 kcat4=0.025;
 KM4=15;
 V5=0.75;
 KM5=15;
 V6=0.75;
 KM6=15;
 kcat7=0.025;
 KM7=15;
 kcat8=0.025;
 KM8=15;
 V9=0.5;
 KM9=15;
 V10=0.5;
 KM10=15;


%run simulation
ODEFUN=@mapkddt;
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Refine', 3);


Tend=8000;
[t,S]=ode45(ODEFUN, [0,Tend], [100 0 300 0 0 300 0 0]);


figure(1)
plot(t, S(:,8), 'k', 'LineWidth',3)
xlabel('Time')
ylabel('Concentration')
title('transient')


%generate dose-reponse.  Larger values of M give greater resolution.

M=10;

input=zeros(1,M);
output=zeros(1,M);


for i=1:M
    V1=0.00+(i/M)*0.4;
    [t,S]=ode45(ODEFUN, [0,Tend], [100 0 300 0 0 300 0 0]);

    input(i)=V1;
    output(i)=S(length(t),8);

end

figure(2)
plot(input, output, 'k', 'Linewidth', 3);
title('dose-response')


end


%dynamics
function dS = mapkddt(t,S)

global V1;
global KI;
global KA;
global KM1;
global V2;
global KM2;
global kcat3;
global KM3;
global kcat4;
global KM4;
global V5;
global KM5;
global V6;
global KM6;
global kcat7;
global KM7;
global kcat8;
global KM8;
global V9;
global KM9;
global V10;
global KM10;

%state variables
MKKK=S(1);
MKKKP=S(2);
MKK=S(3);
MKKP=S(4);
MKKPP=S(5);
MK=S(6);
MKP=S(7);
MKPP=S(8);

%reaction rates
r1=V1*MKKK*(1+KA*MKPP)/((1+(KI*MKPP))*(KM1 + MKKK));
r2=V2*MKKKP/(KM2+MKKKP);
r3=kcat3*MKKKP*MKK/(KM3+MKK);
r4=kcat4*MKKKP*MKKP/(KM4+MKKP);
r5=V5*MKKPP/(KM5+MKKPP);
r6=V6*MKKP/(KM6+MKKP);
r7=kcat7*MKKPP*MK/(KM7+MK);
r8=kcat8*MKKPP*MKP/(KM8+MKP);
r9=V9*MKPP/(KM9+MKPP);
r10=V10*MKP/(KM10+MKP);

%dynamics
dMKKK=r2-r1;
dMKKKP=r1-r2;
dMKK=r6-r3;
dMKKP=r3+r5-r4-r6;
dMKKPP=r4-r5;
dMK=r10-r7;
dMKP=r7+r9-r8-r10;
dMKPP=r8-r9;

dS=[dMKKK dMKKKP dMKK dMKKP dMKKPP dMK dMKP dMKPP]';

    
end


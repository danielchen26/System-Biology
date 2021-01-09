% file lac_operon.m
% model of lac operon in E. coli
% from Santillan et al. (2007) Biophys. J. 92, pp. 3830-3842
% Figure 7.7

function lac_operon
%declare parameters
global delta_M
global delta_Y
global delta_L
global a1
global RToverK1
global K2
global c1
global kL
global KML
global kg
global KMg
global Le
global flag

%assign parameter values
delta_M=0.48;
delta_Y=0.03;
delta_L=0.02;
a1=0.29;
K2=2.92*1e6;
RToverK1=213.2;
c1=18.8;
kL=6*1e4;
KML=680;
kg=3.6*1e3;
KMg=7*1e5;
Le=0;

%initial condition. State is [M Y L]
s0=[0.01  0.1  0];
Tend=2500;
ODEFUN=@lacmodelddt;

%flag=1 for time-series plot (time-varying input profile)
flag=1;

%generate simulation
options=odeset('MaxStep', 0.5, 'Refine', 1);
[t,s] = ode15s(ODEFUN,[0,Tend], s0);
 
 
%generate input profile for plotting
 t_timecourse=[0 499.99 500  999.99  1000   1499.99  1500 1999.99 2000   2500 ];
Le_timecourse=[0 0 50 50     100  100   150 150   0    0];

%generate figure 7.7A
figure(1)
hold on
set(gca,'fontsize',14)
plot(t,s(:,2)/4,'k','Linewidth', 3)
hold off
xlabel('Time (min)')
ylabel('\beta-galactosidase concentration (molecules per cell)')
axis([0 Tend 0 100])
h1=gca;
h2 = axes('Position',get(h1,'Position'));
hold on
plot([0,1],[-1,-1],'k', 'Linewidth', 3)
plot(t_timecourse, Le_timecourse, 'k-.', 'Linewidth', 2)
ylabel('External lactose concentration (\mu M)','fontsize', 14)
legend('\beta-galactosidase (b)', 'external lactose (L_e)', 'fontsize',14)
axis([0 Tend 0 300])
set(h2,'YAxisLocation','right','Color','none','XTickLabel',[],'fontsize', 14)
str1(1) = {'A'};
text(-430,300,str1, 'Fontsize', 40)
 


%generate dose-response 
flag=2;

%dynamics for modified model 
ODEFUN2=@lacnofeedbackmodelddt;

%run simulations to steady state over a range of fixed inputs
N=180;
compYVec=zeros(1,N);
compLeVec=zeros(1,N);
pertYVec=zeros(1,N);
pertLeVec=zeros(1,N);
for i=1:N
Le=120*(i-1)/N;
[t,s] = ode15s(ODEFUN,[0,Tend], s0);
compLeVec(i)=Le;
compYVec(i)=s(length(t),2);
[t,s] = ode15s(ODEFUN2,[0,Tend], s0);
pertLeVec(i)=Le;
pertYVec(i)=s(length(t),2);
end

%generate Figure 7.7B`
figure(2)
set(gca,'fontsize',14)
plot(compLeVec,compYVec/4, 'k-', 'LineWidth', 3)
hold on
plot(pertLeVec,pertYVec/4, 'k--','LineWidth', 3)
xlabel('External lactose concentration (\muM)')
ylabel('\beta-galactosidase concentration (molecules per cell)')
axis([0 100 0 100])
legend('Original model', 'Modified model')
str1(1) = {'B'};
text(-15,100,str1, 'Fontsize', 40)
hold off


end

%dynamics of model
function dS=lacmodelddt(t,s)
                     
global delta_M
global delta_Y
global delta_L
global a1
global RToverK1
global K2
global c1
global kL
global KML
global kg
global KMg
global Le
global flag


%time-varying input profile
if flag==1
Le=0;
if t<500
Le=0.00;
else if t<1000
    Le=50;
    else   if t<1500
            Le=100;
        else if t<2000
            Le=150;
            end
        end
    end
end
end


%declare states
M=s(1);
Y=s(2);
L=s(3);

%dependent states
B=Y/4;
A=L;

%dynamics
dMdt=  a1*1/(1+RToverK1*(K2/(K2+A))^4) - delta_M*M;
dYdt=c1*M - delta_Y*Y;
dLdt=kL*Y*Le/(KML+Le) - 2*kg*B*L/(KMg+L) - delta_L*L;
   
dS=[dMdt dYdt dLdt]';
   
end

%dynamcis for modified model
function dS=lacnofeedbackmodelddt(t,s)
                     
global delta_M
global delta_Y
global delta_L
global a1
global K2
global RToverK1
global c1
global kL
global KML
global kg
global KMg
global Le

%set enzyme level at a fixed value
Enz=40;

%declare states
M=s(1);
Y=s(2);
L=s(3);

%dependent states
B=Y/4;
A=L;

%dynamics
dMdt= a1*1/(1+RToverK1*(K2/(K2+A))^4) - delta_M*M;
dYdt=c1*M - delta_Y*Y;
dLdt=kL*4*Enz*Le/(KML+Le) - kg*Enz*L/(KMg+L) - kg*B*L/(KMg+L) - delta_L*L;
   
dS=[dMdt dYdt dLdt]';
   
end



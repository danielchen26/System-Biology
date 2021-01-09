%file symmetric_network.m
%Model of symmetric network from Figure 4.6. This code generates Figures
%4.7, 4.8, 4.9, and 4.19A

function symmetric_network


%declare model parameters
global k1;
global k2;
global k3;
global k4;
global n1;
global n2;

%declare dynamics and set simulation options
ODEFUN=@symmetricddt;
options=odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Refine', 3);

%Set simulation length
Tend=4;


%Two sets of parametrizations are considered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Asymmetric parametrization    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%declare parameter values
k1=20
k2=20;
k3=5
k4=5
n1=1
n2=4


%Figure 4.7A:
%generate time courses from two separate initial conditions
[t1,S1]=ode45(ODEFUN, [0,Tend], [3,1]);
[t2,S2]=ode45(ODEFUN, [0,Tend], [1,3]);
%generate figure
figure(1)
subplot(2,1,1)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 100 40])
set(gcf,'Position',[0 0 300 400])
set(gca,'fontsize',12)
plot(t1, S1(:,1), 'k', t1, S1(:,2), 'k--', 'LineWidth',3)
axis([0 Tend 0 5])
xlabel('Time')
ylabel('Concentration')
legend('S_1', 'S_2')
str1(1) = {'A'};
text(-.7,5,str1, 'Fontsize', 40)

subplot(2,1,2)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 100 40])
set(gcf,'Position',[0 0 300 400])
set(gca,'fontsize',12)
plot(t2, S2(:,1), 'k', t2, S2(:,2), 'k--', 'LineWidth',3)
axis([0 Tend 0 5])
xlabel('Time')
ylabel('Concentration')
legend('S_1', 'S_2')

    
%figure 4.7B

%generate multiple simulation traces
[t1,S1]=ode45(ODEFUN, [0,Tend], [0,0.5]);
[t2,S2]=ode45(ODEFUN, [0,Tend], [0.5,0]);
[t3,S3]=ode45(ODEFUN, [0,Tend], [1.5,0]);
[t4,S4]=ode45(ODEFUN, [0,Tend], [0,1.5]);
[t5,S5]=ode45(ODEFUN, [0,Tend], [2.5,0]);
[t6,S6]=ode45(ODEFUN, [0,Tend], [3.5,0]);
[t7,S7]=ode45(ODEFUN, [0,Tend], [4.5,0.5]);
[t8,S8]=ode45(ODEFUN, [0,Tend], [0.5,4.5]);
[t9,S9]=ode45(ODEFUN, [0,Tend], [4.5,1.5]);
[t10,S10]=ode45(ODEFUN, [0,Tend], [1.5,4.5]);
[t11,S11]=ode45(ODEFUN, [0,Tend], [4.5,2.5]);
[t12,S12]=ode45(ODEFUN, [0,Tend], [2.5,4.5]);
[t13,S13]=ode45(ODEFUN, [0,Tend], [4.5,3.5]);
[t14,S14]=ode45(ODEFUN, [0,Tend], [3.5,4.5]);

%plot simulation traces
figure(2)
set(gca,'fontsize',14)
hold on
plot(S1(:,1),S1(:,2), 'k', 'LineWidth',1)
plot(S2(:,1),S2(:,2),'k','LineWidth',1)
plot(S3(:,1),S3(:,2),'k','LineWidth',1)
plot(S4(:,1),S4(:,2),'k','LineWidth',1)
plot(S5(:,1),S5(:,2),'k','LineWidth',1)
plot(S6(:,1),S6(:,2),'k','LineWidth',1)
plot(S7(:,1),S7(:,2),'k','LineWidth',1)
plot(S8(:,1),S8(:,2),'k','LineWidth',1)
plot(S9(:,1),S9(:,2),'k','LineWidth',1)
plot(S10(:,1),S10(:,2),'k','LineWidth',1)
plot(S11(:,1),S11(:,2),'k','LineWidth',1)
plot(S12(:,1),S12(:,2),'k','LineWidth',1)
plot(S13(:,1),S13(:,2),'k','LineWidth',1)
plot(S14(:,1),S14(:,2),'k','LineWidth',1)
str1(1) = {'B'};
text(-.5,4.4,str1, 'Fontsize', 40)


%generate nullclines
%s1 nullcline
ns12=linspace(0,5,100);
ns11=k1./(k3*(1.+ns12.^n2));
%s2 nullcline (k2/2 since divided by (1+1^0))
ns21=linspace(0,5,100);
ns22=(k2./(1.+ns21.^n1))/k4;
%plot nullclines
plot(ns11, ns12, 'k--', 'LineWidth',3)
plot(ns21, ns22, 'k--', 'LineWidth',3)

axis([0 4.5 0 4.5])
xlabel('Concentration of S_1')
ylabel('Concentration of S_2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Symmetric parametrization    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%declare parameter values
k1=20
k2=20;
k3=5
k4=5
n1=4
n2=4
    

%figure 4.8A

%generate simualtion curves from two initial conditions
[t1,S1]=ode45(ODEFUN, [0,Tend], [3,1]);
[t2,S2]=ode45(ODEFUN, [0,Tend], [1,3]);

%generate figure
figure(3)
subplot(2,1,1)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 100 40])
set(gcf,'Position',[0 0 300 400])
set(gca,'fontsize',12)
plot(t1, S1(:,1), 'k', t1, S1(:,2), 'k--', 'LineWidth',3)
axis([0 Tend 0 5])
xlabel('Time')
ylabel('Concentration')
legend('S_1', 'S_2')
str1(1) = {'A'};
text(-0.60,5,str1, 'Fontsize', 40)

subplot(2,1,2)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 100 40])
set(gcf,'Position',[0 0 300 400])
set(gca,'fontsize',12)
plot(t2, S2(:,1), 'k', t2, S2(:,2), 'k--', 'LineWidth',3)
axis([0 Tend 0 5])
xlabel('Time')
ylabel('Concentration')
legend('S_1', 'S_2')
    

%figure 4.8B

%generate multiple simulation traces
[t1,S1]=ode45(ODEFUN, [0,Tend], [0,0]);
[t2,S2]=ode45(ODEFUN, [0,Tend], [4.5,4.5]);
[t3,S3]=ode45(ODEFUN, [0,Tend], [0,0.1]);
[t4,S4]=ode45(ODEFUN, [0,Tend], [0.1,0]);
[t5,S5]=ode45(ODEFUN, [0,Tend], [0,0.5]);
[t6,S6]=ode45(ODEFUN, [0,Tend], [0.5,0]);
[t7,S7]=ode45(ODEFUN, [0,Tend], [1.5,0]);
[t8,S8]=ode45(ODEFUN, [0,Tend], [0,1.5]);
[t9,S9]=ode45(ODEFUN, [0,Tend], [4.5,0.5]);
[t10,S10]=ode45(ODEFUN, [0,Tend], [0.5,4.5]);
[t11,S11]=ode45(ODEFUN, [0,Tend], [4.5,1.5]);
[t12,S12]=ode45(ODEFUN, [0,Tend], [1.5,4.5]);
[t13,S13]=ode45(ODEFUN, [0,Tend], [4.5,2.5]);
[t14,S14]=ode45(ODEFUN, [0,Tend], [2.5,4.5]);
[t15,S15]=ode45(ODEFUN, [0,Tend], [4.5,3.5]);
[t16,S16]=ode45(ODEFUN, [0,Tend], [3.5,4.5]);

%plot simulation traces
figure(4)
set(gca,'fontsize',14)
hold on
plot(S1(:,1),S1(:,2), 'k:', 'LineWidth',2)
plot(S2(:,1),S2(:,2),'k:','LineWidth',2)
plot(S3(:,1),S3(:,2),'k','LineWidth',1)
plot(S4(:,1),S4(:,2),'k','LineWidth',1)
plot(S5(:,1),S5(:,2),'k','LineWidth',1)
plot(S6(:,1),S6(:,2),'k','LineWidth',1)
plot(S7(:,1),S7(:,2),'k','LineWidth',1)
plot(S8(:,1),S8(:,2),'k','LineWidth',1)
plot(S9(:,1),S9(:,2),'k','LineWidth',1)
plot(S10(:,1),S10(:,2),'k','LineWidth',1)
plot(S11(:,1),S11(:,2),'k','LineWidth',1)
plot(S12(:,1),S12(:,2),'k','LineWidth',1)
plot(S13(:,1),S13(:,2),'k','LineWidth',1)
plot(S14(:,1),S14(:,2),'k','LineWidth',1)
plot(S15(:,1),S15(:,2),'k','LineWidth',1)
plot(S16(:,1),S16(:,2),'k','LineWidth',1)

%generate nullclines
%s1 nullcline
ns12=linspace(0,5,100);
ns11=k1./(k3*(1.+ns12.^n2));
%s2 nullcline (k2/2 since divided by (1+1^0))
ns21=linspace(0,5,100);
ns22=(k2./(1.+ns21.^n1))/k4;
%plot nullclines
plot(ns11, ns12, 'k--', 'LineWidth',3)
plot(ns21, ns22, 'k--', 'LineWidth',3)
axis([0 4.5 0 4.5])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')
str1(1) = {'B'};
text(-0.55,4.5,str1, 'Fontsize', 40)


%Figure 4.9

%generate multiple simulation traces
[t1,S1]=ode45(ODEFUN, [0,4], [0,0]);
[t2,S2]=ode45(ODEFUN, [0,4], [4,4]);
[t3,S3]=ode45(ODEFUN, [0,4], [0.525,0.5]);
[t4,S4]=ode45(ODEFUN, [0,4], [0.5, 0.525]);
[t5,S5]=ode45(ODEFUN, [0,4], [4,3.9]);
[t6,S6]=ode45(ODEFUN, [0,4], [3.9,4]);
[t7,S7]=ode45(ODEFUN, [0,4], [0.75,0.5]);
[t8,S8]=ode45(ODEFUN, [0,4], [0.5, 0.75]);
[t9,S9]=ode45(ODEFUN, [0,4], [2,1.5]);
[t10,S10]=ode45(ODEFUN, [0,4], [1.5, 2]);

%plot simulation traces
figure(5)
set(gca,'fontsize',14)
hold on
plot(S1(:,1),S1(:,2),'k:','LineWidth',2)
plot(S2(:,1),S2(:,2),'k:','LineWidth',2)
plot(S3(:,1),S3(:,2),'k','LineWidth',1)
plot(S4(:,1),S4(:,2),'k','LineWidth',1)
plot(S5(:,1),S5(:,2),'k','LineWidth',1)
plot(S6(:,1),S6(:,2),'k','LineWidth',1)
plot(S7(:,1),S7(:,2),'k','LineWidth',1)
plot(S8(:,1),S8(:,2),'k','LineWidth',1)
plot(S9(:,1),S9(:,2),'k','LineWidth',1)
plot(S10(:,1),S10(:,2),'k','LineWidth',1)

%generate nullclines
%s1 nullcline
ns12=linspace(0,5,100);
ns11=k1./(k3*(1.+ns12.^n2));
%s2 nullcline (k2/2 since divided by (1+1^0))
ns21=linspace(0,5,100);
ns22=(k2./(1.+ns21.^n1))/k4;
%plot nullclines
plot(ns11, ns12, 'k--', 'LineWidth',2)
plot(ns21, ns22, 'k--', 'LineWidth',2)
axis([0.5 2 0.5 2])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')


%Figure 4.19A

figure(6)
subplot(2,2,1)
set(gca,'fontsize',14)
hold on

k1=10;
%s1 nullcline
ns12=linspace(0,8,200);
ns11=k1./(k3*(1.+ns12.^n2));
%s2 nullcline 
ns21=linspace(0,8,200);
ns22=(k2./(1.+ns21.^n1))/k4;
%setup axes
c=[0.6 0.6 0.6;0 0 0]
set(0,'DefaultAxesColorOrder',c)
plot(ns11, ns12, ns21, ns22,'LineWidth',1)
axis([0 4.2 0 4.2])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')


subplot(2,2,2)
set(gca,'fontsize',14)
hold on
k1=16.2;
%s1 nullcline
ns12=linspace(0,8,200);
ns11=k1./(k3*(1.+ns12.^n2));
%s2 nullcline (k2/2 since divided by (1+1^0))
ns21=linspace(0,8,200);
ns22=(k2./(1.+ns21.^n1))/k4;
%setup axes
c=[0.6 0.6 0.6;0 0 0]
set(0,'DefaultAxesColorOrder',c)
plot(ns11, ns12, ns21, ns22, 'LineWidth',1)
axis([0 4.2 0 4.2])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')


subplot(2,2,3)
set(gca,'fontsize',14)
hold on
k1=20;
%s1 nullcline
ns12=linspace(0,8,200);
ns11=k1./(k3*(1.+ns12.^n2));
%s2 nullcline (k2/2 since divided by (1+1^0))
ns21=linspace(0,8,200);
ns22=(k2./(1.+ns21.^n1))/k4;
%setup axes
c=[0.6 0.6 0.6;0 0 0]
set(0,'DefaultAxesColorOrder',c)
plot(ns11, ns12, ns21, ns22, 'LineWidth',1)
axis([0 4.2 0 4.2])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')


subplot(2,2,4)
set(gca,'fontsize',14)
hold on
k1=35;
%s1 nullcline
ns12=linspace(0,8,200);
ns11=k1./(k3*(1.+ns12.^n2));
%s2 nullcline (k2/2 since divided by (1+1^0))
ns21=linspace(0,8,200);
ns22=(k2./(1.+ns21.^n1))/k4;
%setup axes
c=[0.6 0.6 0.6;0 0 0]
set(0,'DefaultAxesColorOrder',c)
plot(ns11, ns12, ns21, ns22, 'LineWidth',1)
axis([0 7.1 0 4.2])
xlabel('S_1 Concentration')
ylabel('S_2 Concentration')


end

%dynamics for symmetric network model
function dS = symmetricddt(t,S)

global k1 ;
global k2;
global k3;
global k4;
global n1;
global n2;


dS=[k1/(1+S(2)^n2) - k3*S(1); k2/(1+S(1)^n1) - k4*S(2)];

    
end

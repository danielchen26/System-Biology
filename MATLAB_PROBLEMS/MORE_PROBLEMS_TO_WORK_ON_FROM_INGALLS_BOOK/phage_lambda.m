% file phage_lambda.m
% model of phage lambda decision switch
% Figure 7.11


function phage_lambda

%declare parameters
global K1
global K2
global K3
global K4
global a
global b
global delta_r
global delta_c

%assign parameter values
K1=1; 
K2=0.1; 
K3=5; 
K4=0.5; 
delta_r=0.02; 
delta_c=0.02;
a=5; 
b=50;



%set simulation parameters
Tend=6000;
ODEFUN=@phagelambdaddt;
options=odeset('MaxStep', 0.5, 'Refine', 1);


%produce nullclines for Figure 7.11A
figure(1)
hold on
myfun1 = @(r,c) (a + 10*a*K1*(r/2)^2)./(1+K1*(r/2)^2+ K1*K2*(r/2)^3+K3*(c/2)+K3*K4*(c/2)^2)- delta_r*(r)
myfun2 = @(r,c) (b + b*K3*(c/2))./(1+K1*(r/2)^2+ K1*K2*(r/2)^3+K3*(c/2)+K3*K4*(c/2)^2) - delta_c*(c)
a1=ezplot(@(r,c) myfun1(r,c), [-10 350 0 350])
setcurve('color','black','Linewidth', 5) 
set(gca, 'fontsize', 14)
a2=ezplot(@(r,c) myfun2(r,c), [-10 350 0 350])
setcurve('color',[0.5 0.5 0.5],'Linewidth', 5) 
legend('r nullcline', 'c nullcline')
xlabel('cI concentration (nM)')
ylabel('cro concentration (nM)')
title('')
str1(1) = {'A'};
text(-33,250,str1, 'Fontsize', 40)

%generate simulation traces for Figrue 7.11A
[t1,s1] = ode15s(ODEFUN,[0,Tend], [250,50]);
[t2,s2] = ode15s(ODEFUN,[0,Tend], [250,150]);
[t3,s3] = ode15s(ODEFUN,[0,Tend], [250,250]);
[t4,s4] = ode15s(ODEFUN,[0,Tend], [150,250]);
[t5,s5] = ode15s(ODEFUN,[0,Tend], [100,250]);
[t6,s6] = ode15s(ODEFUN,[0,Tend], [50,250]);
[t7,s7] = ode15s(ODEFUN,[0,Tend], [25,250]);
[t8,s8] = ode15s(ODEFUN,[0,Tend], [17,250]);
[t9,s9] = ode15s(ODEFUN,[0,Tend], [3,2]);
[t10,s10] = ode15s(ODEFUN,[0,Tend], [6,25]);

plot(s1(:,1), s1(:,2),'--', 'Linewidth', 2)
plot(s2(:,1), s2(:,2),'--','Linewidth', 2)
plot(s3(:,1), s3(:,2),'--','Linewidth', 2)
plot(s4(:,1), s4(:,2),'--','Linewidth', 2)
plot(s5(:,1), s5(:,2),'--','Linewidth', 2)
plot(s6(:,1), s6(:,2),'--','Linewidth', 2)
plot(s7(:,1), s7(:,2),'--','Linewidth', 2)
plot(s8(:,1), s8(:,2),'--','Linewidth', 2)
plot(s9(:,1), s9(:,2),'--','Linewidth', 2)
plot(s10(:,1), s10(:,2),'--','Linewidth', 2)

axis([0 250 0 250])

hold off

%generate nullclines for Figure 7.11B
delta_cI=0.2;
figure(2)
hold on
myfun1 = @(r,c) (a + 10*a*K1*(r/2)^2)./(1+K1*(r/2)^2+ K1*K2*(r/2)^3+K3*(c/2)+K3*K4*(c/2)^2)- delta_r*(r)
myfun2 = @(r,c) (b + b*K3*(c/2))./(1+K1*(r/2)^2+ K1*K2*(r/2)^3+K3*(c/2)+K3*K4*(c/2)^2) - delta_c*(c)
a1=ezplot(@(r,c) myfun1(r,c), [-10 350 0 350])
setcurve('color','black','Linewidth', 5) 
set(gca, 'fontsize', 14)
a2=ezplot(@(r,c) myfun2(r,c), [-10 350 0 350])
setcurve('color',[0.5 0.5 0.5],'Linewidth', 5) 
legend('r nullcline', 'c nullcline')
xlabel('cI concentration (nM)')
ylabel('cro concentration (nM)')
title('')
str1(1) = {'B'};
text(-33,250,str1, 'Fontsize', 40)

%generate simulation traces for Figrue 7.11A
[t1,s1] = ode15s(ODEFUN,[0,50], [250,15]);
[t2,s2] = ode15s(ODEFUN,[0,Tend], [250,65]);
[t3,s3] = ode15s(ODEFUN,[0,Tend], [250,100]);
[t4,s4] = ode15s(ODEFUN,[0,Tend], [250,150]);
[t5,s5] = ode15s(ODEFUN,[0,Tend], [250,200]);
[t6,s6] = ode15s(ODEFUN,[0,Tend], [250,250]);
[t7,s7] = ode15s(ODEFUN,[0,Tend], [50,250]);
[t8,s8] = ode15s(ODEFUN,[0,Tend], [150,250]);
[t9,s9] = ode15s(ODEFUN,[0,Tend], [15,20]);
[t10,s10] = ode15s(ODEFUN,[0,Tend], [15,5]);

plot(s1(:,1), s1(:,2),'--', 'Linewidth', 2)
plot(s2(:,1), s2(:,2),'--','Linewidth', 2)
plot(s3(:,1), s3(:,2),'--','Linewidth', 2)
plot(s4(:,1), s4(:,2),'--','Linewidth', 2)
plot(s5(:,1), s5(:,2),'--','Linewidth', 2)
plot(s6(:,1), s6(:,2),'--','Linewidth', 2)
plot(s7(:,1), s7(:,2),'--','Linewidth', 2)
plot(s8(:,1), s8(:,2),'--','Linewidth', 2)
plot(s9(:,1), s9(:,2),'--','Linewidth', 2)
plot(s10(:,1), s10(:,2),'--','Linewidth', 2)

axis([0 250 0 250])
hold off

end

%dynamics
function dS=phagelambdaddt(t,s)
                     
global K1
global K2
global K3
global K4
global a
global b
global delta_r
global delta_c

%declare states
r=s(1);
c=s(2);
rd=r/2;
cd=c/2;

%dynamics
drdt= (a + 10*a*K1*rd^2)/(1+K1*rd^2+ K1*K2*rd^3+K3*cd+K3*K4*cd^2)- delta_r*r;
dcdt= (b + b*K3*cd)/(1+K1*rd^2+ K1*K2*rd^3+K3*cd+K3*K4*cd^2) - delta_c*c;
      
dS=[drdt dcdt]';
   
end

%change properties of last curve in current figure
%Examples:
%     setcurve('color','red')
%     setcurve('color','green','linestyle','--')
%Type  help plot  to see available colors and line styles 
function setcurve(varargin)
h=get(gca,'children');
set(h(1),varargin{:})
end

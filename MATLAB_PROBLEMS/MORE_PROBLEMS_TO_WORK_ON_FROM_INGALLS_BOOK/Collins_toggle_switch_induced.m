% file Collins_toggle_switch_induced.m
% model of Collins toggle switch
% from Gardiner et al. (2000) Nature 403, pp. 339-342
% Figures 7.14B

function toggleswitch


%declare parameters
global a1
global a2
global beta
global gamma
global i1
global i2

%assign parameter values
a1=3;
a2=2.5;
beta=4;
gamma=4;


%set simulation parameters
x0=[0.075    2.5]';
Tend=50
ODEFUN=@toggleswitchddt;
options=odeset('Refine', 6);


%Simulations for Figure 7.14B
i1=0;
i2=10;
    
[t1,s1] = ode15s(ODEFUN,[0,Tend], [0,0]);
[t2,s2] = ode15s(ODEFUN,[0,Tend], [1,0]);
[t3,s3] = ode15s(ODEFUN,[0,Tend], [2,0]);
[t4,s4] = ode15s(ODEFUN,[0,Tend], [3,0]);
[t5,s5] = ode15s(ODEFUN,[0,Tend], [3.5,0]);
[t6,s6] = ode15s(ODEFUN,[0,Tend], [3.5,0.5]);
[t7,s7] = ode15s(ODEFUN,[0,Tend], [3.5,1.5]);
[t8,s8] = ode15s(ODEFUN,[0,Tend], [3.5,2.5]);
[t9,s9] = ode15s(ODEFUN,[0,Tend], [3.5,3.5]);
[t10,s10] = ode15s(ODEFUN,[0,Tend], [3,3.5]);
[t11,s11] = ode15s(ODEFUN,[0,Tend], [2,3.5]);
[t12,s12] = ode15s(ODEFUN,[0,Tend], [1,3.5]);
[t13,s13] = ode15s(ODEFUN,[0,Tend], [0,3]);
[t14,s14] = ode15s(ODEFUN,[0,Tend], [0,2]);
[t15,s15] = ode15s(ODEFUN,[0,Tend], [0,1]);

figure(1)
hold on

%generate nullclines
myfun1 = @(p1,p2) a1/(1+(p2/(1+i2))^beta) - p1
myfun2 = @(p1,p2) a2/(1+(p1/(1+i1))^gamma) - p2
a1=ezplot(@(p1,p2) myfun1(p1,p2), [0 3.5 0 3.5])
setcurve('color','black','Linewidth', 5) 
set(gca, 'fontsize', 14)
a2=ezplot(@(p1,p2) myfun2(p1,p2), [0 3.5 0 3.5])
setcurve('color',[0.5 0.5 0.5],'Linewidth', 5) 
legend('repressor 1 nullcline', 'repressor 2 nullcline')
xlabel('repressor 1 concentration (arbitrary units)')
ylabel('repressor 2 concentration (arbitrary units)')
title('')
str1(1) = {'B'};
text(-0.4,3.5,str1, 'Fontsize', 40)

%plot simulation traces
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
plot(s11(:,1), s11(:,2),'--','Linewidth', 2)
plot(s12(:,1), s12(:,2),'--','Linewidth', 2)
plot(s13(:,1), s13(:,2),'--','Linewidth', 2)
plot(s14(:,1), s14(:,2),'--','Linewidth', 2)
plot(s15(:,1), s15(:,2),'--','Linewidth', 2)

axis([0 3.5 0 3.5])

hold off





end

%dynamics
function dS=toggleswitchddt(t,x)


global a1
global a2
global beta
global gamma
global i1
global i2


            

dS =[a1/(1+(x(2)/(1+i2))^beta) - x(1),  a2/(1+(x(1)/(1+i1))^gamma) - x(2)]';

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

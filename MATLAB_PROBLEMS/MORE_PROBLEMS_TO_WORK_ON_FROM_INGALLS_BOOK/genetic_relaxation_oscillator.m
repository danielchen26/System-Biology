%file genetic_relaxation_oscillator.m
%Hasty synthetic oscillator model
%from Hasty et al. (2002) Phys. Rev. Lett. 88 148101
%Figure 7.23
%NOTE: typos in the published book: 
%page 263, caption of figure 7.23, the value for gamma_y should be 0.012, not 0.06.
%page 303 Problem 7.8.14, the value of gamma_y should be  0.15, (not 1 as written).


function relaxation_oscillator

clear all

%declare model parameters
global alpha
global sigma
global gammax
global gammay
global ay

%assign parameter values
alpha=11;
sigma=2;
gammax=0.2;
gammay=0.012;
ay=0.2;

%set simulation parameters
ODEFUN=@hastyddt;
options=odeset('Refine', 6);
Tend=300;


%set initial condition
x0=[0.3963    2.3346]';


%run simulation
[t,s]=ode45(ODEFUN, [0 Tend], x0, options);

%produce figure 7.23A
 figure(1)
 set(gca,'fontsize',14)
 plot(t,s(:,1),'k', t,s(:,2), 'k--','Linewidth', 3)
 axis([0 300 0 3])
 xlabel('Time (arbitrary units)')
 ylabel('Concentration (arbitrary units)')
 legend('activator X', 'repressor Y')
 str1(1) = {'A'};
 text(-30,3,str1, 'Fontsize', 40)
 

%produce figure 7.23B
figure(2)
 
%generate nullclines
hold on
myfun1 = @(x,y) (1+x^2+alpha*sigma*x^4)/((1+x^2+sigma*x^4)*(1+y^4)) - gammax*x
myfun2 = @(x,y) (ay)*((1+x^2+alpha*sigma*x^4)/((1+x^2+sigma*x^4)*(1+y^4))) - gammay*y
a1=ezplot(@(x,y) myfun1(x,y), [0 1.5 0 3])
setcurve('color','black','Linewidth', 5) 
set(gca, 'fontsize', 14)
a2=ezplot(@(x,y) myfun2(x,y), [0 1.5 0 3])
setcurve('color',[0.5 0.5 0.5],'Linewidth', 5) 
legend('x nullcline (activator)', 'y nullcline (repressor)')
xlabel('activator concentration (x)')
ylabel('repressor concentration (y)')
title('')
str1(1) = {'B'};
text(-.2,3,str1, 'Fontsize', 40)


%generate simulation along limit cycle
 x0=[0.3963    2.3346]';
times=linspace(0, 74, 100);
[t,s]=ode45(ODEFUN, times, x0, options);
plot(s(:,1), s(:,2), 'k.', 'Linewidth', 3)
axis([0 1.5 1.5 3])

end


%model dynamics
function dS=hastyddt(t,s)

global alpha
global sigma
global gammax
global gammay
global ay


x=s(1);
y=s(2);

xddt= (1+x^2+alpha*sigma*x^4)/((1+x^2+sigma*x^4)*(1+y^4)) - gammax*x;
yddt= ay*((1+x^2+alpha*sigma*x^4)/((1+x^2+sigma*x^4)*(1+y^4))) - gammay*y;
            
dS =[xddt,  yddt]';

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


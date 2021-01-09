% file Collins_toggle_switch.m
% model of Collins toggle switch
% from Gardiner et al. (2000) Nature 403, pp. 339-342
% Figures 1.7, 7.13, 7.14, 7.15

function toggleswitch


%declare parameters
global a1
global a2
global beta
global gamma
global transient_flag
global i1
global i2

%assign parameter values
a1=3;
a2=2.5;
beta=4;
gamma=4;
%no inputs
i1=0;
i2=0;


%set simulation parameters
x0=[0.075    2.5]';
Tend=50
ODEFUN=@toggleswitchddt;
options=odeset('Refine', 6);

%Generate Figures 1.7 and 7.13
transient_flag=1; %this sets a time-varying input profile
[t,s]=ode45(ODEFUN, [0,Tend], x0, options);
figure(1)
set(gca,'fontsize',14)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 100 40])
set(gcf,'Position',[0 0 1000 400])
plot(t,s(:,1), 'k', t,s(:,2), 'k--', 'Linewidth', 3)
axis([0 Tend 0 3.5])
xlabel('Time')
ylabel('Concentration')
leg=legend('Protein 1', 'Protein 2');
str1(1)={'inducer 2 active'}
text(11,3.25,str1, 'Fontsize', 14)
str2(1)={'inducer 1 active'}
text(31,3.25,str2, 'Fontsize', 14)



 
% generate phaseportraits

transient_flag=0;



%Simulations for Figure 7.14A

   
[t1,s1] = ode15s(ODEFUN,[0,Tend], [0,0]);
[t2,s2] = ode15s(ODEFUN,[0,Tend], [0.5,0]);
[t3,s3] = ode15s(ODEFUN,[0,Tend], [1,0]);
[t4,s4] = ode15s(ODEFUN,[0,Tend], [1.5,0]);
[t5,s5] = ode15s(ODEFUN,[0,Tend], [2,0]);
[t6,s6] = ode15s(ODEFUN,[0,Tend], [2.5,0]);
[t7,s7] = ode15s(ODEFUN,[0,Tend], [3,0]);
[t8,s8] = ode15s(ODEFUN,[0,Tend], [3.5,0]);
[t9,s9] = ode15s(ODEFUN,[0,Tend], [3.5,0.5]);
[t10,s10] = ode15s(ODEFUN,[0,Tend], [3.5,1]);
[t11,s11] = ode15s(ODEFUN,[0,Tend], [3.5,1.5]);
[t12,s12] = ode15s(ODEFUN,[0,Tend], [3.5,2]);
[t13,s13] = ode15s(ODEFUN,[0,Tend], [3.5,2.5]);
[t14,s14] = ode15s(ODEFUN,[0,Tend], [3.5,3]);
[t15,s15] = ode15s(ODEFUN,[0,Tend], [3.5,3.5]);
[t16,s16] = ode15s(ODEFUN,[0,Tend], [3,3.5]);
[t17,s17] = ode15s(ODEFUN,[0,Tend], [2.5,3.5]);
[t18,s18] = ode15s(ODEFUN,[0,Tend], [2,3.5]);
[t19,s19] = ode15s(ODEFUN,[0,Tend], [1.5,3.5]);
[t20,s20] = ode15s(ODEFUN,[0,Tend], [1,3.5]);
[t21,s21] = ode15s(ODEFUN,[0,Tend], [0.5,3.5]);
[t22,s22] = ode15s(ODEFUN,[0,Tend], [0,3.5]);
[t23,s23] = ode15s(ODEFUN,[0,Tend], [0,3]);
[t24,s24] = ode15s(ODEFUN,[0,Tend], [0,2.5]);
[t25,s25] = ode15s(ODEFUN,[0,Tend], [0,2]);
[t26,s26] = ode15s(ODEFUN,[0,Tend], [0,1.5]);
[t27,s27] = ode15s(ODEFUN,[0,Tend], [0,1]);
[t28,s28] = ode15s(ODEFUN,[0,Tend], [0,0.5]);

figure(2)
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
str1(1) = {'A'};
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
plot(s16(:,1), s16(:,2),'--','Linewidth', 2)
plot(s17(:,1), s17(:,2),'--','Linewidth', 2)
plot(s18(:,1), s18(:,2),'--','Linewidth', 2)
plot(s19(:,1), s19(:,2),'--','Linewidth', 2)
plot(s20(:,1), s20(:,2),'--','Linewidth', 2)
plot(s21(:,1), s21(:,2),'--','Linewidth', 2)
plot(s22(:,1), s22(:,2),'--','Linewidth', 2)
plot(s23(:,1), s23(:,2),'--','Linewidth', 2)
plot(s24(:,1), s24(:,2),'--','Linewidth', 2)
plot(s25(:,1), s25(:,2),'--','Linewidth', 2)
plot(s26(:,1), s26(:,2),'--','Linewidth', 2)
plot(s27(:,1), s27(:,2),'--','Linewidth', 2)
plot(s28(:,1), s28(:,2),'--','Linewidth', 2)


axis([0 3.5 0 3.5])

hold off


%Generate Figure 7.15
%This is a "brute force" parameter sweep

%set mesh size.  N=100 in the figure in the text, but that takes a long
%time to run.

N=8; %100


%Three cases are considered: for the three values of beta, gamma
%%%%Case 1
beta=4;
gamma=4;

grid=zeros(N);  
tolgrid=grid;

for i=1:N
    for j=1:N
     
%sweep logarithmically through a1-a2 space.
a1=exp(0+6*i/N);
a2=exp(0+6*j/N);

%Run simulations from low p1/high p2 and low p2/high p1
[t1,s1] = ode15s(ODEFUN,[0,2000], [100,0]);
[t2,s2] = ode15s(ODEFUN,[0,2000], [0,100]);

%record the difference in s1 at steady state
grid(i,j)=abs(s1(length(t1),1)-s2(length(t2),1));


if grid(i,j)<0.001 %if the two simulations come to the same point
tolgrid(i,j)=0;   %tolgrid is set to zero if the difference is negligible
else
    tolgrid(i,j)=grid(i,j); %otherwise, the difference is recorded.
end

    end
end


%For each x-coordinate (i.e. a2), determine the y-coordinate (i.e. a1) at which bistability
%begins
xcurve1=zeros(N,2);

for xi=1:N
    xcurve1(xi,1)=xi;
    yi=1;
    while (tolgrid(xi,yi)==0)&&(yi<N)
        yi=yi+1;
    end
    xcurve1(xi,2)=yi;
end

%For each y-coordinate (i.e. a1), determine the x-coordinate (i.e. a2) at which bistability
%begins
ycurve1=zeros(N,2);

for yi=1:N
    ycurve1(yi,1)=yi;
    xi=1;
    while (tolgrid(xi,yi)==0)&&(xi<N)
        xi=xi+1;
    end
    ycurve1(yi,2)=xi;
end

%plot the resulting points in Figure 7.15
figure(4)
hold on
plot(xcurve1(:,1), xcurve1(:,2), 'k.')
plot(ycurve1(:,2), ycurve1(:,1), 'k.')


%%%%Case 2
beta=2;
gamma=2;
grid=zeros(N);
tolgrid=grid;

for i=1:N
    for j=1:N
        
a1=exp(0+6*i/N);
a2=exp(0+6*j/N);

[t1,s1] = ode15s(ODEFUN,[0,2000], [100,0]);
[t2,s2] = ode15s(ODEFUN,[0,2000], [0,100]);

grid(i,j)=abs(s1(length(t1),1)-s2(length(t2),1));

if grid(i,j)<0.001
tolgrid(i,j)=0;
else
    tolgrid(i,j)=grid(i,j);
end

    end
end


xcurve2=zeros(N,2);

for xi=1:N
    xcurve2(xi,1)=xi;
    yi=1;
    while (tolgrid(xi,yi)==0)&&(yi<N)
        yi=yi+1;
    end
    xcurve2(xi,2)=yi;
end

i=1;
while xcurve2(i,2)==N
    i=i+1;
end
xcurve2=xcurve2(i:N,:);

ycurve2=zeros(N,2);

for yi=1:N
    ycurve2(yi,1)=yi;
    xi=1;
    while (tolgrid(xi,yi)==0)&&(xi<N)
        xi=xi+1;
    end
    ycurve2(yi,2)=xi;
end

i=1;
while ycurve2(i,2)==N
    i=i+1;
end
ycurve2=ycurve2(i:N,:);

hold on
plot(xcurve2(:,1), xcurve2(:,2), 'kx')
plot(ycurve2(:,2), ycurve2(:,1), 'kx')

%%%%Case 3
beta=1.5;
gamma=1.5;
grid=zeros(N);
tolgrid=grid;

for i=1:N
    for j=1:N
        
a1=exp(0+6*i/N);
a2=exp(0+6*j/N);

[t1,s1] = ode15s(ODEFUN,[0,2000], [100,0]);
[t2,s2] = ode15s(ODEFUN,[0,2000], [0,100]);

grid(i,j)=abs(s1(length(t1),1)-s2(length(t2),1));

if grid(i,j)<0.001
tolgrid(i,j)=0;
else
    tolgrid(i,j)=grid(i,j);
end

    end
end

xcurve3=zeros(N,2);

for xi=1:N
    xcurve3(xi,1)=xi;
    yi=1;
    while (tolgrid(xi,yi)==0)&&(yi<N)
        yi=yi+1;
    end
    xcurve3(xi,2)=yi;
end

i=1;
while xcurve3(i,2)==N
    i=i+1;
end
xcurve3=xcurve3(i:N,:);

ycurve3=zeros(N,2);

for yi=1:N
    ycurve3(yi,1)=yi;
    xi=1;
    while (tolgrid(xi,yi)==0)&&(xi<N)
        xi=xi+1;
    end
    ycurve3(yi,2)=xi;
end

i=1;
while ycurve3(i,2)==N
    i=i+1;
end
ycurve3=ycurve3(i:N,:);
    
    
%figure(4)

plot(xcurve3(:,1), xcurve3(:,2), 'k*')
plot(ycurve3(:,2), ycurve3(:,1), 'k*')

set(gca, 'Fontsize', 14)
xlabel('log(\alpha_2)')
ylabel('log(\alpha_1)')
%To reset tickmarks:
 set(gca,'YTickLabel','0|1|2|3|4|5|6|7')
 set(gca,'XTickLabel','0|1|2|3|4|5|6|7')


end

%dynamics
function dS=toggleswitchddt(t,x)


global a1
global a2
global beta
global gamma
global transient_flag
global i1
global i2

%time-varying input profile for Figure 7.13
if transient_flag>0
if t< 10
    i1=0;i2=0;
else if t<20
        i1=0;i2=10;
else if t<30
       i1=0;i2=0;
    else if t<40
            i1=10;i2=0;
        else if t<=50
               i1=0;i2=0;
            end
        end
    end
    end
end
end
            

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


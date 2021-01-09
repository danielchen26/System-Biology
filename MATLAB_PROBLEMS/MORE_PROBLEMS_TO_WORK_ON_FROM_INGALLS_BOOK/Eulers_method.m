%File eulers_method.m
%Figure 2.7 Numerical Simulation by Euler's Method

%Final time:
tf=2;

%generate exact solution a(t)
t=linspace(0,tf,100);
a=exp(-t);

%Solution b(t) via Euler's method with step-size h=2/3
%step-size:
h=2/3;
%initial condition:
b(1)=1;
%number of steps:
N=tf/h;
%mesh:
t2=linspace(0,tf,N+1);
%Euler's method:
for i=2:N+1
b(i)=b(i-1) - h*b(i-1);
end

%Solution c(t) via Euler's method with step-size h=1/3
%step-size:
h=1/3;
%initial condition:
c(1)=1;
%number of steps:
N=tf/h;
%mesh:
t3=linspace(0,tf,N+1);
%Euler's method:
for i=2:N+1
c(i)=c(i-1) - h*c(i-1);
end

%generate plot
figure(1)
set(gca,'fontsize',14)
hold on
plot(t,a, 'k', 'Linewidth', 3)
plot(t2,b, 'k--', t3, c, 'k:', 'Linewidth', 2)
plot(t2,b, 'ko' ,t3, c, 'ks', 'MarkerSize',15)
axis([0 tf 0 1])
xlabel('Time (arbitrary units)')
ylabel('Concentration (arbitrary units)')


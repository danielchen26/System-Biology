%file stochastic_decay.m
%stochastic implementation (Gillespie SSA) of a collection of spontaneously
%decaying molecules.
%Figure 7.38

%the network is simply  X -> . with propensity a=k N_x, where N_x is the
%number of molecules of X

%declare vectors for storing the steps in molecule count and time
X=zeros(1,1000);
t=zeros(1,1000);

%set parameter value
k=1;

%set initial condition for molecule count and for simulation time.
X(1)=10;
t(1)=0;

%run SSA steps until 5 time units have passed
i=1;
while t < 5

%calculate propensity
if X(i)>0
a=k*X(i);
end

%update counter
i=i+1;

%update time
t(i)=t(i-1)+log(1/rand(1))/a;

%decrement state
X(i)=max(X(i-1)-1,0);


end


X=X(1:i);
t=t(1:i);

%produce Figure 7.38A
figure(1)
subplot(1,3,3)
stairs(t,X, 'k', 'linewidth', 2)
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')
axis([0 5 0 10])
str1(1) = {'C'};
hold on
%%add deterministic plot
det_t=linspace(0,5);
sol=10*exp(-det_t);
plot(det_t, sol, 'k--')

%repeat for Figure 7.38B

X=zeros(1,1000);
t=zeros(1,1000);
X(1)=100;
t(1)=0;
i=1;

while t < 5

%calculate propensity
if X(i)>0
a=k*X(i);
end

%update counter
i=i+1;

%update time
t(i)=t(i-1)+log(1/rand(1))/a;

%decrement state
X(i)=max(X(i-1)-1,0);


end


X=X(1:i);
t=t(1:i);

figure(1)
subplot(1,3,2)
stairs(t,X, 'k', 'linewidth', 2)
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')
axis([0 5 0 100])
str1(1) = {'B'};
hold on
%%add deterministic plot
det_t=linspace(0,5);
sol=100*exp(-det_t);
plot(det_t, sol, 'k--')

%repeat for Figure 7.38C

X=zeros(1,1000);
t=zeros(1,1000);
X(1)=1000;
t(1)=0;
i=1;

while t < 5

%calculate propensity
if X(i)>0
a=k*X(i);
end

%update counter
i=i+1;

%update time
t(i)=t(i-1)+log(1/rand(1))/a;

%decrement state
X(i)=max(X(i-1)-1,0);


end


X=X(1:i);
t=t(1:i);

hold on
figure(1)
subplot(1,3,1)
stairs(t,X, 'k', 'linewidth', 2)
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')
axis([0 5 0 1000])
str1(1) = {'A'};
hold on
%add deterministic plot
det_t=linspace(0,5);
sol=1000*exp(-det_t);
plot(det_t, sol, 'k--')
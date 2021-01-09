%file stochastic_decay_ensemble.m
%stochastic implementation (Gillespie SSA) of a collection of spontaneously
%decaying molecules.
%Figure 7.45

%the network is simply  X -> . with propensity a=k N_x, where N_x is the
%number of molecules of X

clear all

%set parameter value
k=1;

%set ensemble size
N=3;

%set initial condition for molecule count and for simulation time, for each member of ensemble
X1=zeros(10,N);
t1=zeros(10,N);
X1(1,:)=10;
t1(1,:)=0;

%generate ensemble
for j=1:N

%run ten steps of the simulation algorithm
for i=1:10
    
%calculate propensity
if X1(i,j)>0
a=k*X1(i,j);
end

%update counter
i=i+1;

%update time
t1(i,j)=t1(i-1,j)+log(1/rand(1))/a;

%decrement state
X1(i,j)=max(X1(i-1,j)-1,0);


end


end

%collect vector of reaction time events 
rxntimes=[];
for j=1:N
rxntimes=union(rxntimes, t1(:,j));
end

%interpolate molecule counts for plotting (otherwise the ensemble of staircase graphs lie
%on top of one another)
for j=1:N
X1interp(:,j)=interp1(t1(:,j), X1(:,j), rxntimes);
end

X1interp(isnan(X1interp)) = 0;
X1mean=mean(X1interp');

%produce Figure 7.45A
figure(1)
hold on

subplot(1,3,1)
colormap(gray)
c=[0 0 0;0.4 0.4 0.4;0.8 0.8 0.8]
set(0,'DefaultAxesColorOrder',c)
hold on

plot(rxntimes, X1mean, 'k', 'Linewidth', 3)
for j=1:N
plot(t1(:,j), X1(:,j), 'color', [0.5 0.5 0.5])
end
set(gca,'fontsize',12)

axis([0 5 0 10])
ylabel('Number of Molecules')
xlabel('Time')

%repeat for larger ensemble sizes

%set ensemble size
N=10
X2=zeros(10,N);
t2=zeros(10,N);
X2(1,:)=10;
t2(1,:)=0;

for j=1:N

for i=1:10
    
%calculate propensity
if X2(i,j)>0
a=k*X2(i,j);
end

%update counter
i=i+1;

%update time
t2(i,j)=t2(i-1,j)+log(1/rand(1))/a;

%decrement state
X2(i,j)=max(X2(i-1,j)-1,0);


end


end

rxntimes=[];
for j=1:N
rxntimes=union(rxntimes, t2(:,j));
end


for j=1:N
X2interp(:,j)=interp1(t2(:,j), X2(:,j), rxntimes);
end

X2interp(isnan(X2interp)) = 0;
X2mean=mean(X2interp');

%produce Figure 7.45B
figure(1)
hold on
subplot(1,3,2)
hold on
for j=1:N
plot(t2(:,j), X2(:,j), 'color', [0.5 0.5 0.5])
end
plot(rxntimes, X2mean, 'k','Linewidth', 3)
axis([0 5 0 10])
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')


%set ensemble size
N=30
X3=zeros(10,N);
t3=zeros(10,N);
X3(1,:)=10;
t3(1,:)=0;

for j=1:N

for i=1:10
    
%calculate propensity
if X3(i,j)>0
a=k*X3(i,j);
end

%update counter
i=i+1;

%update time
t3(i,j)=t3(i-1,j)+log(1/rand(1))/a;

%decrement state
X3(i,j)=max(X3(i-1,j)-1,0);


end


end

rxntimes=[];
for j=1:N
rxntimes=union(rxntimes, t3(:,j));
end


for j=1:N
X3interp(:,j)=interp1(t3(:,j), X3(:,j), rxntimes);
end

X3interp(isnan(X3interp)) = 0;
X3mean=mean(X3interp');

%produce Figure 7.45C
figure(1)
hold on
subplot(1,3,3)
hold on
for j=1:N
plot(t3(:,j), X3(:,j), 'color', [0.5 0.5 0.5])
end
plot(rxntimes, X3mean, 'k','Linewidth', 3)
axis([0 5 0 10])
set(gca,'fontsize',12)
ylabel('Number of Molecules')
xlabel('Time')


%File gillespie.m
%stochastic simulation of -> s1 -> s2 ->

clear all

%set parameter value
v0=20;
k1=1;
k2=2;

%set number of reaction
Tend=10000;


%set initial condition for concentrations, and for time t 

X=zeros(Tend, 2); t=zeros(Tend,1)
X(1,:)=[0 0]; t(1)=0;


for j=1:Tend

%calculate propensities
a1=v0;
a2=k1*X(j,1);
a3=k2*X(j,2);

asum=a1+a2+a3;


%select nect reaction
mu=rand(1);
z1=0;z2=0;z3=0;
if 0 <= mu && mu < a1/asum
z1=1;
    else if a1/asum <= mu && mu  < (a1+a2)/asum
        z2=1;
        else 
            z3=1;
        end
end

%update state

X(j+1,1)= X(j,1)+z1-z2;
X(j+1,2)=X(j,2)+z2-z3;
           
    
%update time
t(j+1)=t(j)+log(1/rand(1))/asum;
   
%update counter
j=j+1;

end



figure(1)

stairs(t, X(:,1), 'k', 'linewidth', 2)
hold on
stairs(t, X(:,2), 'g', 'linewidth', 2)
ylabel('Number of Molecules')
xlabel('Time')



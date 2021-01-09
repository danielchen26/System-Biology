%file stochastic_gene_expression.m
%stochastic implementation (Gillespie SSA) of a constitutive gene
%expression
%Figure 7.46


clear all

%set kinetic parameter values
kr=10;
kp=6;
gr=1;
gp=1;


%set initial condition for molecule counts (M, P) and for simulation time

X=[0 0 0];

%set number of reaction
Tend=10000;

for j=1:Tend

%calculate propensities
a1=kr;
a2=kp*X(j,2);
a3=gr*X(j,2);
a4=gp*X(j,3);

asum=a1+a2+a3+a4;

%update counter
j=j+1;

%update time
X(j,1)=X(j-1,1)+log(1/rand(1))/asum;

%state transition
mu=rand(1);
if 0 <= mu && mu < a1/asum
X(j,2)=X(j-1,2)+1;
X(j,3)=X(j-1,3);
    else if a1/asum <= mu && mu  < (a1+a2)/asum
        X(j,2)=X(j-1,2);
        X(j,3)=X(j-1,3)+1;
        else if (a1+a2)/asum <= mu && mu < (a1+a2+a3)/asum
            X(j,2) = max(X(j-1,2)-1,0);
            X(j,3)=X(j-1,3);
            else
            X(j,2)=X(j-1,2);
            X(j,3) = max(X(j-1,3)-1,0);
            end
        end
     end
end



figure(1)
set(gca,'fontsize',14)

%produce figure 7.46A
subplot(1,2,1)
set(gca,'fontsize',14)
stairs(X(:,1), X(:,2), 'k', 'linewidth', 2)
hold on
stairs(X(:,1), X(:,3), 'color', [0.5 0.5 0.5], 'linewidth', 2)
hold on
axis([0 70 0 140])
ylabel('Number of Molecules')
xlabel('Time')
legend('mRNA', 'protein')
str1(1) = {'A'};
text(-18,140,str1, 'Fontsize', 40)
hold off


%modify burst size
b=5;

kr=10/b;
kp=6;
gr=1;
gp=1;


%set initial condition for concentrations, and for time t, 
X=[0 0 0];

%set number of reaction
Tend=10000;

for j=1:Tend

%calculate propensities
a1=kr;
a2=kp*X(j,2);
a3=gr*X(j,2);
a4=gp*X(j,3);

asum=a1+a2+a3+a4;



%update counter
j=j+1;

%update time
X(j,1)=X(j-1,1)+log(1/rand(1))/asum;

%state transition
mu=rand(1);
if 0 <= mu && mu < a1/asum
X(j,2)=X(j-1,2)+b;
X(j,3)=X(j-1,3);
    else if a1/asum <= mu && mu  < (a1+a2)/asum
        X(j,2)=X(j-1,2);
        X(j,3)=X(j-1,3)+1;
        else if (a1+a2)/asum <= mu && mu < (a1+a2+a3)/asum
            X(j,2) = max(X(j-1,2)-1,0);
            X(j,3)=X(j-1,3);
            else
            X(j,2)=X(j-1,2);
            X(j,3) = max(X(j-1,3)-1,0);
            end
        end
     end
end


%produce Figure 7.46B
figure(1)
subplot(1,2,2)
stairs(X(:,1), X(:,2), 'k', 'linewidth', 2)
hold on
stairs(X(:,1), X(:,3), 'color', [0.5 0.5 0.5], 'linewidth', 2)
hold on
axis([0 70 0 140])
set(gca,'fontsize',14)
ylabel('Number of Molecules')
xlabel('Time')
legend('mRNA', 'protein')
str1(1) = {'B'};
text(-18,140,str1, 'Fontsize', 40)
hold off


% File chain1.m  
% Model of two-species reaction chain (-> s1 -> s2 ->)}
% with irreversible mass action kinetics}  
%
% The first function generates the simulation and produces a plot:
function chain1
% 
% Declaration of model parameters.  
% Global declarations allow parameters to be shared between the two functions 
global v0; 
global k1; 
global k2; 
%
% Declaration of model equations.  The function chain1dtt is defined in the second part of this file.
%
ODEFUN=@chain1ddt;
%
% Assignment of parameter values 
%
v0=5; 
k1=3; 
k2=1; 
%
% Set initial conditions as a vector 
s0=[1,0];  
%
% Specify the length of time to be simulated  
%
Tend=5;    % The simulation will run over the interval [0, Tend]} 
%
% Generate the simulation.  
%
[t,s]= ode45(ODEFUN, [0,Tend], s0);
%The vector t now contains the time-points on which the simulation is defined.
%The matrix s has a row for each time-point in t.  Each row contains the vector [s1,s2]  at that time-point
%
% Open a plotting window  
%
figure(1);
%
% Plot the time-course for s1 (with a thicker-than default curve) 
%
plot(t, s(:,1), 'LineWidth',3)
%
% Retain the curve in the window  while others are plotted  
%
hold on; 
%
% Plot the time-courses for s2 as a dashed curve
%
plot(t, s(:,2), '--', 'LineWidth',3)
%
% Specify the window size ([0, Tend]x[0,5]) 
%
axis([0 Tend 0 5]);
%
% Add axis labels and a legend 
%
xlabel('Time'); 
ylabel('Concentration'); 
legend('s\_1', 's\_2');
%
%
% Open a second plotting window  
%
figure(2); 
%
% Plot the time-course for the rate of the second reaction 
%
plot(t, k1*s(:,1), 'LineWidth',3)
%
% Specify the window size ([0, 2]x[0,5]) 
%
axis([0 2 0 7]);
%
% Add axis labels and a legend
%
xlabel('Time'); 
ylabel('Rate');  
legend('rate of reaction two'); 
%
% Close the function chain1:
%
end 
% 
% 
%  
% 
% The second function defines the model equations.  This function has input (t,s) -- the time and state at a given point, and output ds -- the right-hand-side of the system of differential equations  
function ds=chain1ddt(t,s)
% 
% These global declarations allow these parameters to be shared between the two functions
global v0; 
global k1; 
global k2; 
% 
% Assignment of state variables 
% 
s1=s(1); 
s2=s(2);
% 
% Declaration of model equations
% 
ds1 = v0 - k1*s1; 
ds2 = k1*s1 - k2*s2;
% 
% Assignment of the output vector 
% 
ds = [ds1;ds2];
% 
% Close the function chain1ddt:
%
end  
% 

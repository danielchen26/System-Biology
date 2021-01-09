%file Hodgkin_Huxley.m
%Hodgkin-Huxley model of excitable barnacle muscle fiber
%reviewed in Rinzel (1990) Bulletin of Mathematical Biology 52 pp. 5-23.
%Figure 1.9 and problem 8.6.4

function Hodgkin_Huxley

%declare model parameters
global E_N
global E_K
global g_N_bar
global g_K_bar
global g_leak
global E_leak
global C_M

%assign parameter values
E_N=55;
E_K=-72;
g_N_bar = 120.0;
g_K_bar = 36.0;
g_leak = 0.30;
E_leak = -49.0;   
C_M = 1.0;


%initial condition. state = [V m h n]
s0=[-59.8977    0.0536    0.5925    0.3192];


%run simulation
Tend=100;
ODEFUN=@HHmodelddt;
options=odeset('MaxStep', 0.05, 'Refine', 3);
[t,s] = ode15s(ODEFUN,[0,Tend], s0, options);
 
%produce figure 1.9
figure(1)
set(gca,'fontsize',14)
plot(t,s(:,1), 'k', 'Linewidth', 2)
axis([0 Tend -80 50])
xlabel('Time (msec)')
ylabel('Membrane voltage (mV)')



end

%dynamics
function dS=HHmodelddt(t,S)
  

global E_N
global E_K
global g_N_bar
global g_K_bar
global g_leak
global E_leak
global C_M

V=S(1);
m=S(2);
h=S(3);
n=S(4)
 
m_inf = (0.10*(V+35)/(1.0 - exp(-(V+35)/10.0)))/((0.10*(V+35)/(1.0 - exp(-(V+35)/10.0)))+(4.0*exp( -(V+60)/18.0)));
t_m   = 1.0/((0.10*(V+35)/(1.0 - exp(-(V+35)/10.0)))+(4.0*exp( -(V+60)/18.0)));
dmdt  = (m_inf - m)/t_m;
   
h_inf = (0.07*exp(-(V+60)/20))/((0.07*exp(-(V+60)/20))+(1/(exp(-(V+30)/10)+1)));
t_h   = 1.0/((0.07*exp(-(V+60)/20))+(1/(exp(-(V+30)/10)+1)));
dhdt  = (h_inf - h)/t_h;

I_N = g_N_bar*(V-E_N)*m^3*h;  

t_n = 1.0/((0.01*(V+50)/(1.0 - exp(-(V+50)/10.0)))+(0.125*exp( -(V+60)/80)));
n_inf = (0.01*(V+50)/(1.0 - exp(-(V+50)/10.0)))/((0.01*(V+50)/(1.0 - exp(-(V+50)/10.0)))+(0.125*exp( -(V+60)/80)));
dndt  = (n_inf - n)/t_n;

I_K = g_K_bar*(V-E_K)*n^4; 

I_leak = g_leak*(V-E_leak);
   
  
   Istim=0;
   if t>20
       Istim=  -6.65;
   end
   if t>21
       Istim=0;
   end
   if t>60
       Istim=-6.85;
   end
   if t>61
       Istim=0;
   end
   
   dVdt = (- I_N - I_K - I_leak - Istim)/C_M;
   
   
   dS=[dVdt dmdt dhdt dndt]';
   
end

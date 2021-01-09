function  fit_lotka_volterra

%we take y1(t) to be prey population and y2(t) to be predator population
%vector of parameters is (p(1), p(2), p(3), p(4)) = (a, b, c, r).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LVPERUN: MATLAB script M-file to run Lotka-Volterra
%parameter estimation example.

% minimizes the function lverr with
% MATLAB’s built-in nonlinear minimizing routine fminsearch
% using a good guess from a quick search using another method

global H
global L
global years

lvdata
years=0:1:20;

% copied from curve_fitting_matlab_howard_2009_tamu_m469.pdf
% in above notes, guess = [.47; .024; .023; .76];
% found using linear regression on approximatio of logarithmic derivatives
%figure lv047_0024_0023_076.*
%initial_guess = 0.4700    0.0240    0.0230    0.7600
%found_parameters = 0.5486    0.0283    0.0264    0.8375
%error_for_these = 744.7935

guess = [1; .02; .02; 1];
%figure lv1_002_002_1.*
%initial_guess =  1.0000    0.0200    0.0200    1.0000
%found_parameters =  0.5484    0.0283    0.0265    0.8384
%error_for_these =  744.8950

%guess = [0.3; 0.3; 0.3 ; 0.3];
%figure lv1_02_03_03_03.*
%initial_guess =  0.3000    0.3000    0.3000    0.3000
%found_parameters =  0.8808    0.0743    0.0087    0.0569
%error_for_these =  2.2942e+04

[p,error]=fminsearch(@lverr, guess);

initial_guess = guess'

found_parameters = p'

error_for_these = error

% error is sum of squares

[t,y]=ode45(@lvpe,years,[H(1);L(1)],[],p);
subplot(2,1,1)
plot(t,y(:,1),years,H,'o')
subplot(2,1,2)
plot(t,y(:,2),years,L,'o')

% homework: try other initial guesses!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = lvpe(t,y,p)

%LVPE: ODE for example Lotka-Volterra parameter
%estimation example. p(1)=a, p(2) = b, p(3) = c, p(4) = r.

value=[p(1)*y(1)-p(2)*y(1)*y(2);-p(4)*y(2)+p(3)*y(1)*y(2)];
end %lvpe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function error = lverr(p)

%LVERR: Function defining error function for
%example with Lotka-Volterra equations.

% take as input valuesthe parameter vector p and return squared error E(p)
% ode statement passes to MATLAB a vector of times at which we have observations
%years=[0 , 20]

[t,y] = ode45(@lvpe,years,[H(1);L(1)],[],p);

value = (y(:,1)-H').^2+(y(:,2)-L').^2;
% primes transpose data vectors H and L
error = sum(value);
end %lverr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lvdata
  H=[30 47.2 70.2 77.4 36.3 20.6 18.1 21.4 22 25.4 27.1 ...
     40.3 57 76.6 52.3 19.5 11.2 7.6 14.6 16.2 24.7];
  L=[4 6.1 9.8 35.2 59.4 41.7 19 13 8.3 9.1 7.4 ...
     8 12.3 19.5 45.7 51.1 29.7 15.8 9.7 10.1 8.6];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %program

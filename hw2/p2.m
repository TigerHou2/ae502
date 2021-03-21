%% AE 502 HW2 Problem 2, Spring 2021
%   Tiger Hou
close all
clear;clc

mu = 3.986e14; % m^3/s^2, Earth gravitational parameter
J2 = 1082.63e-6; % J2 perturbation coefficient
req = 6378.137e3; % km, equatorial radius

% initial conditions
a = 7000e3;
e = 0.05;
i = deg2rad(45);
o = deg2rad(0);
w = deg2rad(45);
M0 = deg2rad(0);
E0 = kepler(M0,e);
f0 = 2 * atan(sqrt((1+e)/(1-e))*tan(E0/2));

% define the perturbation equation in the inertial frame as p(r)
p = @(r) -3/2 * J2 * (mu/norm(r)^2) * (req/norm(r))^2 * ...
    [ ( 1-5*(r(3)/norm(r))^2 ) * r(1)/norm(r); ...
      ( 1-5*(r(3)/norm(r))^2 ) * r(2)/norm(r); ...
      ( 3-5*(r(3)/norm(r))^2 ) * r(3)/norm(r) ];

% calculate orbit period
T = 2*pi * sqrt(a^3/mu); % seconds

% ode45
param0 = [a,e,i,o,w,M0]';
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[tOut,paramOut] = ode45(@(t,params)ff(params,mu,p),[0,10*T],param0,options);

% plotting
plot(tOut,paramOut(:,5))

%% function definitions
function param_dot = ff(params,mu,A)
% takes the 6x1 orbit parameters and computes the derivative using Gauss'
% variation of parameters
%   the orbit parameters are a, e, i, o, w, M
%   also applies perturbation in the form of a function handle A
%   which takes argument A(params)

a = params(1);
e = params(2);
i = params(3);
o = params(4);
w = params(5);
M = params(6);
E = kepler(M,e);
f = 2 * atan(sqrt((1+e)/(1-e))*tan(E/2));
p = a*(1-e^2); % semi-latus rectum
n = sqrt(mu/a^3); % mean motion
b = sqrt(a*p);

[R,V] = Get_Orb_Vects([a,e,i,o,w,f],mu);
r = norm(R);
H = cross(R,V);
h = norm(H);

% find the LVLH reference frame basis vectors
ir = R / r;
in = H / h;
it = cross(in,ir);

aR = A(R)' * ir; % radial perturbation
aT = A(R)' * it; % theta perturbation
aN = A(R)' * in; % normal perturbation

dadt = (2*a^2/h) * ( e*sin(f)*aR + p/r*aT );
dedt = 1/h * ( p*sin(f)*aR + ((p+r)*cos(f)+r*e)*aT );
didt = 1/h * r*cos(w+f) * aN;
dodt = (r*sin(w+f)) / (h*sin(i)) * aN;
dwdt = 1/h/e * ( -p*cos(f)*aR + (p+r)*sin(f)*aT ) ...
        - (r*sin(w+f)*cos(i)) / (h*sin(i)) * aN;
dmdt = n + b/(a*h*e) * ( (p*cos(f)-2*r*e)*aR - (p+r)*sin(f)*aT );

param_dot = [dadt;dedt;didt;dodt;dwdt;dmdt];

end

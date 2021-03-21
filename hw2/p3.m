%% AE 502 HW2 Problem 3, Spring 2021
%   Tiger Hou
close all
clear;clc
latexify

%% Setup
% Earth orbit general parameters
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

% find initial position and velocity vectorsr0
[r0,v0] = Get_Orb_Vects([a,e,i,o,w,f0],mu);

% define the J2 perturbation equation as p_J2(r)
p_J2 = @(r) -3/2 * J2 * (mu/norm(r)^2) * (req/norm(r))^2 * ...
    [ ( 1-5*(r(3)/norm(r))^2 ) * r(1)/norm(r); ...
      ( 1-5*(r(3)/norm(r))^2 ) * r(2)/norm(r); ...
      ( 3-5*(r(3)/norm(r))^2 ) * r(3)/norm(r) ];

% Vallado exponential drag model (taking only relevant portions)
drag_Vallado = [200e3,  250e3, 2.789e-10, 37.105e3; ...
                250e3,  300e3, 7.248e-11, 45.546e3; ...
                300e3,  350e3, 2.418e-11, 53.628e3; ...
                350e3,  400e3, 9.518e-12, 53.298e3; ...
                400e3,  450e3, 3.725e-12, 58.515e3; ...
                450e3,  500e3, 1.585e-12, 60.828e3; ...
                500e3,  600e3, 6.967e-13, 63.822e3; ...
                600e3,  700e3, 1.454e-13, 71.835e3; ...
                700e3,  800e3, 3.614e-14, 88.667e3; ...
                800e3,  900e3, 1.170e-14, 124.64e3; ...
                900e3, 1000e3, 5.245e-15, 181.05e3];
rho = @(r) sum(...
           ( (r-req) >= drag_Vallado(:,1) ) ...
        .* ( (r-req) <  drag_Vallado(:,2) ) ...
        .* drag_Vallado(:,3) ...
        .* exp(-(r-req-drag_Vallado(:,1)) ./ drag_Vallado(:,4)) ...
               ) ...
        /  sum(... this line checks if at least one drag model is matched
           ( (r-req) >= drag_Vallado(:,1) ) ... otherwise division by zero
        .* ( (r-req) <  drag_Vallado(:,2) ) );
% drag model
Cd = 2.0;
A = 5; % m^2
m = 600; % kg
% define the drag perturbation equation as p_drag(rv)
p_drag = @(rv) -1/2 * Cd * A / m * rho(norm(rv(1:3))) ...
            * norm(rv(4:6)) * rv(4:6);

% define overall perturbation model
p = @(rv) p_J2(rv(1:3)) + p_drag(rv);

% calculate orbit period
T = 2*pi * sqrt(a^3/mu); % seconds

%% propagate for J2 + drag case
% ode45
rv0 = [r0;v0];
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[tDrag,rvDrag] = ode45(@(t,rv)ff(rv,mu,p),[0,10*T],rv0,options);

% plotting
rvDrag = rvDrag';
plot3(rvDrag(1,:), rvDrag(2,:),rvDrag(3,:))
axis equal

% convert position, velcoity data into orbit parameters
N = size(rvDrag,2);
paramDrag = nan(size(rvDrag));
for j = 1:N
    [a_,e_,i_,o_,w_,f_] = Get_Orb_Params(rvDrag(1:3,j),rvDrag(4:6,j),mu);
    E_ = 2 * atan(sqrt((1-norm(e_))/(1+norm(e_)))*tan(f_/2));
    M_ = E_ - norm(e_)*sin(E_);
    paramDrag(:,j) = [a_,norm(e_),i_,o_,w_,M_]';
end

%% propagate for J2 (without drag)
% ode45
rv0 = [r0;v0];
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[tNoDrag,rvNoDrag] = ode45(@(t,rv)ff(rv,mu,p_J2),tDrag,rv0,options);

% plotting
rvNoDrag = rvNoDrag';
plot3(rvNoDrag(1,:), rvNoDrag(2,:),rvNoDrag(3,:))
axis equal

% convert position, velcoity data into orbit parameters
N = size(rvNoDrag,2);
paramNoDrag = nan(size(rvNoDrag));
for j = 1:N
    [a_,e_,i_,o_,w_,f_] = Get_Orb_Params(rvNoDrag(1:3,j),rvNoDrag(4:6,j),mu);
    E_ = 2 * atan(sqrt((1-norm(e_))/(1+norm(e_)))*tan(f_/2));
    M_ = E_ - norm(e_)*sin(E_);
    paramNoDrag(:,j) = [a_,norm(e_),i_,o_,w_,M_]';
end

%% Compare results
ylabel_vec = {'Semi-Major Axis, m', ...
              'Eccentricity', ...
              'Inclination, deg', ...
              'RAAN, deg', ...
              'Argument of Periapsis, deg', ...
              'Mean Anomaly, deg'};
xlabel_val = 'Time, s';
for j = 1:6
    subplot(3,2,j)
    delta = paramDrag(j,:)-paramNoDrag(j,:);
    if j == 6 % correction for mean anomaly loopback from 2*pi to 0
        delta = mod(delta,2*pi);
    end
    if j >= 3 % conversion to degrees
        delta = rad2deg(delta);
    end
    plot(tDrag,delta)
    xlabel(xlabel_val)
    ylabel(ylabel_vec{j})
    setgrid
end
latexify(20,20)

%% function definitions
function rv_dot = ff(rv,mu,p)
% takes the 6x1 position & velocity vector and computes the derivative
%   also applies perturbation in the form of a function handle p
%   which takes argument p(rv)

rv_dot = nan(6,1);

% velocity
rv_dot(1:3) = rv(4:6);
% acceleration
rv_dot(4:6) = -mu/norm(rv(1:3))^3*rv(1:3) + p(rv);

end

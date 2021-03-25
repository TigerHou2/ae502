%% AE 502 HW2 Problem 1, Spring 2021
%   Tiger Hou
close all
clear;clc

%% Part a

tic

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

% define the perturbation equation as p(rv)
p = @(rv) -3/2 * J2 * (mu/norm(rv)^2) * (req/norm(rv))^2 * ...
    [ ( 1-5*(rv(3)/norm(rv))^2 ) * rv(1)/norm(rv); ...
      ( 1-5*(rv(3)/norm(rv))^2 ) * rv(2)/norm(rv); ...
      ( 3-5*(rv(3)/norm(rv))^2 ) * rv(3)/norm(rv) ];
  
% calculate orbit period
T = 2*pi * sqrt(a^3/mu); % seconds

% ode45
rv0 = [r0;v0];
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[tOut,rvOut] = ode45(@(t,rv)ff(rv,mu,1,p),[0,10*T],rv0,options);

toc

% plotting
rvOut = rvOut';
figure(1)
plot3(rvOut(1,:), rvOut(2,:),rvOut(3,:),'r','LineWidth',1.2)
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')
axis equal
setgrid
latexify(16,14,14)

%% Part b
% convert position, velcoity data into orbit parameters
N = size(rvOut,2);
paramOut = nan(size(rvOut));
for j = 1:N
    [a_,e_,i_,o_,w_,f_] = Get_Orb_Params(rvOut(1:3,j),rvOut(4:6,j),mu);
    E_ = 2 * atan(sqrt((1-norm(e_))/(1+norm(e_)))*tan(f_/2));
    M_ = E_ - norm(e_)*sin(E_);
    paramOut(:,j) = [a_,norm(e_),i_,o_,w_,M_]';
end

% plot results
ylabel_vec = {'SMA (a), m', ...
              'ECC (e)', ...
              'INC (i), deg', ...
              'RAAN ($\Omega$), deg', ...
              'AOP ($\omega$), deg', ...
              'MA (M), deg'};
xlabel_val = 'Orbit Count';
figure(2)
for j = 1:6
    subplot(3,2,j)
    dat = paramOut(j,:);
    if j == 4 % RAAN loop-around after 2*pi
        dat = mod(dat+pi,2*pi)-pi;
    end
    if j == 6 % mean anomaly loop-around after 2*pi
        dat = mod(dat,2*pi);
    end
    if j >= 3 % conversion to degrees
        dat = rad2deg(dat);
    end
    plot(tOut/T,dat,'Linewidth',1.2)
    xlabel(xlabel_val)
    ylabel(ylabel_vec{j})
    setgrid
end
latexify(20,18,14)

%% Part c
n = sqrt(mu/a^3); % mean motion
p = a*(1-e^2); % semi-latus rectum
dadt = 0;
dedt = 0;
didt = 0;
dodt = -3/2*J2*n*(req/p)^2*cos(i);
dwdt = 3/4*J2*n*(req/p)^2*(5*cos(i)^2-1);
dM0dt = 3/4*J2*n*(req/p)^2*sqrt(1-e^2)*(3*cos(i)^2-1);
dMdt = dM0dt + n;

param0 = [a,e,i,o,w,M0]';

paramMean = ([dadt,dedt,didt,dodt,dwdt,dMdt]') * (tOut') + param0;

for j = 1:6
    subplot(3,2,j)
    dat = paramOut(j,:)-paramMean(j,:);
    if j == 4 % RAAN loop-around after 2*pi
        dat = mod(dat+pi,2*pi)-pi;
    end
    if j == 6 % mean anomaly loop-around after 2*pi
        dat = mod(dat+pi,2*pi)-pi;
    end
    if j >= 3 % conversion to degrees
        dat = rad2deg(dat);
    end
%     hold on
    plot(tOut/T,dat,'Linewidth',1.2)
    setgrid
    grid minor
%     hold off
    xlabel(xlabel_val)
    ylabel(ylabel_vec{j})
    if j == 4 % add legend in the relatively empty plot
%         legend('Numerical','Mean Motion','Location','best')
    end
end
latexify(20,18,14)

% figure(2)
% plot(tOut,paramMean(6,:))
% hold on
% plot(tOut,paramOut(6,:))
% hold off

%% function definitions
function rv_dot = ff(rv,mu,N,p)
% takes the 6xN position & velocity vector and computes the derivative
%   where N is the number of particles to track
%   also applies perturbation in the form of a function handle p
%   which takes argument p(r)

rv_dot = zeros(6,N);

for i = 1:N
    % velocity
    rv_dot(1:3,i) = rv(4:6,i);
    % acceleration
    rv_dot(4:6,i) = -mu/norm(rv(1:3,i))^3*rv(1:3,i) + p(rv(1:3,i));
end

end

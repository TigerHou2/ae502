%% AE 502 HW3 Problem 2, Spring 2021
%   Tiger Hou

close all
clear;clc
latexify

%% problem setup

% Earth orbit general parameters
mu = 3.986e14; % m^3/s^2, Earth gravitational parameter
J2 = 1082.63e-6; % J2 perturbation coefficient
req = 6378.137e3; % m, equatorial radius

% orbit parameters
alt_geo = 35786e3; % m, GEO altitude
a = alt_geo + req;
e = 0;
i = 0;
o = 0;
w = 0;
f = 0;

m0 = 400; % kg, initial mass of spacecraft
thrust = 0.136; % N, low thrust
Isp = 3100; % s, specific impulse
g0 = 9.81; % m/s^2, Earth gravity
c = Isp * g0; % exhaust velocity

% find initial position and velocity vectors
[r0,v0] = Get_Orb_Vects([a,e,i,o,w,f],mu);

% define the low thrust engine as a perturbation
% with current mass as the fourth state
pT = @(v,m) [ v / norm(v) * thrust / m; -thrust/c ];

% define the perturbation equation as p(r,v,m,t)
p = @(r,v,m) -3/2 * J2 * (mu/norm(r)^2) * (req/norm(r))^2 * ...
    [ ( 1-5*(r(3)/norm(r))^2 ) * r(1)/norm(r); ...
      ( 1-5*(r(3)/norm(r))^2 ) * r(2)/norm(r); ...
      ( 3-5*(r(3)/norm(r))^2 ) * r(3)/norm(r); ...
      0 ] + ...
    pT(v,m);
  
% calculate orbit period
T = 2*pi * sqrt(a^3/mu); % seconds


%% propagate orbit

% ode45
rv0 = [r0;v0;m0];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tOut,rvOut] = ode45(   @(t,rv)ff(rv,mu,p),...
                        [0,24*3600*100],rv0,options);
rvOut = rvOut';
                    
%% check escape
                    
% check escape velocity
escaped = floor(vecnorm(rvOut(4:6,:)) ...
                ./ sqrt( 2*mu ./ vecnorm(rvOut(1:3,:))));
            
idx = find(escaped,1,'first');
if ~isempty(idx)
    disp([  'The spacecraft escaped at t = ' ...
            num2str(tOut(idx)/24/3600) ' days.'])
else
    disp([  'The spacecraft did not escape after' ...
            num2str(tOut(end)/24/3600) ' days.'])
end

%% plot results

% plotting
figure(1)
plot3(rvOut(1,:), rvOut(2,:),rvOut(3,:),'r','LineWidth',1.2)
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')
% axis equal
setgrid
latexify(16,14,14)


%% plot orbit parameters

% convert position, velcoity data into orbit parameters
N = idx-1;
paramOut = nan(size(rvOut,1),N);
for j = 1:N
    [a_,e_,i_,o_,w_,f_] = Get_Orb_Params(rvOut(1:3,j),rvOut(4:6,j),mu);
    E_ = 2 * atan(sqrt((1-norm(e_))/(1+norm(e_)))*tan(f_/2));
    M_ = E_ - norm(e_)*sin(E_);
    paramOut(:,j) = [a_,norm(e_),i_,o_,w_,M_,rvOut(7,j)]';
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
    plot(tOut(1:N)/T,dat,'Linewidth',1.2)
    xlabel(xlabel_val)
    ylabel(ylabel_vec{j})
    setgrid
end
latexify(20,18,14)


%% function definitions
function rv_dot = ff(rv,mu,p)
% takes the 6xN position & velocity vector and computes the derivative
%   where N is the number of particles to track
%   also applies perturbation in the form of a function handle p
%   which takes argument p(r)

rv_dot = zeros(7,1);

% velocity
rv_dot(1:3) = rv(4:6);
% acceleration and mass
rv_dot(4:7) = - mu/norm(rv(1:3))^3*[rv(1:3);0] ...
                + p(rv(1:3), rv(4:6), rv(7));

end
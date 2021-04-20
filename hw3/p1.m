%% AE 502 HW3 Problem 1, Spring 2021
%   Tiger Hou

close all
clear;clc
latexify

%% Part a

% Mars general parameters
mu = 4.282837e4; % km^3/s^2, Mars gravitational parameter
J2 = 1960.45e-6; % J2 perturbation coefficient
req = 3389.5; % km, equatorial radius
Omars = 7.0879e-5; % rad/s, Mars spin rate
P = 687 * 24 * 3600; % seconds, Mars orbit period around Sun

% orbit parameters
e = 0.1;
alt = 1000; % km, orbit altitude
a = alt + req;

n = sqrt(mu/a^3); % mean motion
p = a*(1-e^2); % semi-latus rectum

dodt_sso = 2*pi / P; % orbit precession rate of SSO at Mars
i_sso = acos( dodt_sso / (-3/2*J2*n*(req/p)^2) );

disp(['The sun-synchronous inclination is ' ...
        num2str(rad2deg(i_sso)) ' deg.'])

%% Part b

% initial conditions
i = i_sso;
o = deg2rad(0);
w = deg2rad(0);
M0 = deg2rad(0);
E0 = kepler(M0,e);
f0 = 2 * atan(sqrt((1+e)/(1-e))*tan(E0/2));

% find initial position and velocity vectors
[r0,v0] = Get_Orb_Vects([a,e,i,o,w,f0],mu);

% define the perturbation equation as p(rv)
p = @(rv) -3/2 * J2 * (mu/norm(rv)^2) * (req/norm(rv))^2 * ...
    [ ( 1-5*(rv(3)/norm(rv))^2 ) * rv(1)/norm(rv); ...
      ( 1-5*(rv(3)/norm(rv))^2 ) * rv(2)/norm(rv); ...
      ( 3-5*(rv(3)/norm(rv))^2 ) * rv(3)/norm(rv) ];
  
% calculate orbit period
T = 2*pi * sqrt(a^3/mu); % seconds

% final time
tf = 10 * 24 * 3600;

% ode45
rv0 = [r0;v0];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[tOut,rvOut] = ode45(   @(t,rv)ff(rv,mu,1,p,t,Omars),...
                        [0,tf],rv0,options);

% plotting
rvOut = rvOut';
figure(1)
colorplot(rvOut(1,:), rvOut(2,:),rvOut(3,:),...
            'colormap','bone',...
            'linewidth',1.2)
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')
axis equal
setgrid
latexify(16,14,14)

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
latexify(26,16,14)


%% Part c

v_PCPF = nan(size(rvOut(4:6,:)));

for j = 1:size(v_PCPF,2)
%     theta = tOut(j) * Omars;
%     R_eci2ecef = [  cos(theta), sin(theta), 0; ...
%                    -sin(theta), cos(theta), 0; ...
%                     0         , 0         , 1];
%     v_PCPF(:,j) = R_eci2ecef * rvOut(4:6,j);
    v_PCPF(:,j) = rvOut(4:6,j) + cross(Omars*[0;0;1],rvOut(1:3,j));
end

rr = vecnorm(rvOut(1:3,:));
p2 = rvOut(3,:) ./ rr;

E = 1/2 * vecnorm(v_PCPF).^2 ...
  - 1/2 * Omars^2 * sum(rvOut(1:2,:).^2) ...
  - mu ./ rr ...
       .* ( 1 - (req./rr).^2 .* J2 .* (3*p2.^2-1)/2 );

E0 = E(1);

figure;
loglog(tOut, abs( (E-E0)/E0 ), 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('$\frac{||E(t)-E(t_0)||}{||E(t_0)||}$')
setgrid
latexify(18,14,16)


%% function definitions
function rv_dot = ff(rv,mu,N,p,t,Omars)
% takes the 6xN position & velocity vector and computes the derivative
%   where N is the number of particles to track
%   also applies perturbation in the form of a function handle p
%   which takes argument p(r)

rv_dot = zeros(6,N);

for i = 1:N
    theta = t*Omars;
    R_eci2ecef = [  cos(theta), sin(theta), 0; ...
                   -sin(theta), cos(theta), 0; ...
                    0         , 0         , 1];
    R_ecef2eci = R_eci2ecef';
    % velocity
    rv_dot(1:3,i) = rv(4:6,i);
    % acceleration
    rv_dot(4:6,i) = -mu/norm(rv(1:3,i))^3*rv(1:3,i) ...
                  + R_ecef2eci*p(R_eci2ecef*rv(1:3,i));
end

end

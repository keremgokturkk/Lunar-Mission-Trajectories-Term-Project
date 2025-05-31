% LUNAR MISSION TRAJECTORIES | TERM PROJECT | K.G. ULUTAS
% PART 4 - ASCENT
clear, clc, close all;
global m G n mu_moon Rm T_max mt
%CONSTANTS
n=2;               % Number of the objects;
G = 6.6742867e-11; % Universal Gravitational Constant [N*m*m/kg/kg]
Re=6.3781366e6;  % Earth's radius
m2=5.9721424e24; %  Earth's mass
m1=7.3457576e22; %  Moon's mass
m3=1000;         %  Spacecraft's mass
mp = 6352.2081;
mt = m3 + mp;
Rm = 1737400;
mu_earth=G*m2;
mu_moon=G*m1;
v_park_moon = sqrt(mu_moon/2097400); % Corrected parking orbit velocity
Isp = 336; % Monomethyl Hydrazin
g0 = 9.80665;
T_max = 15000; % Maximum thrust in Newtons
% t=0:0.01:2*pi; %
%m1=2e20; %kg
%m2=5e20; %kg
L=2.0974e6; %m
m=[m1 mt]; % Mass vector
Lm1=mt/(m1+mt)*L;
Lm2=m1/(mt+mt)*L;
w=sqrt(G*(m1+mt)/L^3);
b01=[-Lm1 0 0];        % Initial position for m1
v01=[0 -w*Lm1 0];      % Initial velocity for m1
b02=[-1.236948993031341e+06 -1.220047386657381e+06 0]; % Initial position for m2 which is the surface of Moon
v02=[0 0 0];       % Increased initial velocity for m2
% Initial Conditions Vector
b0=[b01 v01 b02 v02];
% ÇÖZÜM ARALIĞI
tspan1=10; % Initial liftoff for a short time
t_initial=[0 tspan1]; % [s]
options = odeset('RelTol',1e-9,'AbsTol', 1e-9,'Events',@liftoff_start);
[t1,y1,te,ye,ie] = ode45(@SystemsOfEquations,t_initial,b0,options);
if ~isempty(te)
    b0_current = ye;
    T = t1;
    Y = y1;
    
    % ------------------- thrust loop patch -----------------------
a_quad   = 1.53785e-3;             % s^-2, from your derivation
tau      = 1;                  % s, time constant for radial speed damping

for c = 1:1200
    mt = mt - 4.5523045;
    % 1) extract state
    r_vec     = b0_current(1,7:9)' ;   % [m]
    v_vec     = b0_current(1,10:12)';   % [m/s]
    r_mag     = norm(r_vec);
    nu        = -r_vec/r_mag;          % outward unit vector

    % 2) actual & reference radial speeds (+outward)
    v_r       = dot(v_vec, nu);
    t_now     = c;                   % seconds since liftoff
    v_r_ref   = (a_quad*t_now^2 - 1200*a_quad*t_now);

    % 3) desired net radial accel to drive error → 0
    a_net_des = -(v_r - v_r_ref)/tau;  % +outward if v_r<v_r_ref

    % 4) lunar gravity (inward positive)
    a_grav    = mu_moon/r_mag^2;       % m/s²

    % 5) engine radial accel required
    a_eng_rad = a_net_des + a_grav;    % m/s² (+outward)

    % 6) saturate to thruster capability
    a_max     =  T_max/mt;
    a_eng_rad = max(-a_max, min(a_eng_rad, a_max));

    % 7) radial thrust vector (N)
    thrust_radial = a_eng_rad * mt * nu;  % outward if +, else inward

    % 8) unchanged tangential thrust
    h_vec      = cross(r_vec, v_vec);
    z_hat      = h_vec/norm(h_vec);
    t_hat      = cross(z_hat, r_vec);
    t_hat      = t_hat/norm(t_hat);
    thrust_tangential = t_hat * T_max*0.3684869;

    % 9) apply both accelerations
    thrust_accel_rad = thrust_radial    / mt;
    thrust_accel_tan = thrust_tangential / mt;
    b0_current(1,10:12) = b0_current(1,10:12) ...
                        + thrust_accel_rad' ...
                        + thrust_accel_tan';

    % 10) propagate 1 s
    [t_step,y_step] = ode45(@SystemsOfEquations, [0 1], b0_current(1,:), options);
    b0_current      = y_step(end,:);
    t_step          = t_step + T(end);
    T = [T; t_step];
    Y = [Y; y_step];
end
% -------------------------------------------------------------
period_circ2 = 2*pi*sqrt((2097400^3) / mu_moon);
tspan_f = T(end)+period_circ2;
[tstep_f,ystep_f] = ode45(@SystemsOfEquations, [1 period_circ2], b0_current(1,:), options);

T = [T; tstep_f];
Y = [Y; ystep_f];
norm(v_vec)
% norm(r_vec);

else
    fprintf('Liftoff start event did not occur.\n');
end
position = [Y(end,7) Y(end,8)];
velocity = [Y(end,10) Y(end,11)];
Vector = [position, velocity];
%DO NOT CHANGE ANY CODE FROM THE PLOT SECTION
figure
circle = 0:pi/5000:2*pi;
hold on
plot(Y(:,1),Y(:,2),'r',...
Y(:,7),Y(:,8),'b',...
Y(length(T),1)+Rm*cos(circle),Y(length(T),2)+Rm*sin(circle),'r',...
Y(length(T),1)+(Rm+360000)*cos(circle),Y(length(T),2)+(Rm+360000)*sin(circle),'r--');
% Get the final position and velocity vectors.
final_position = [Y(end,7), Y(end,8)];
final_velocity = [Y(end,10), Y(end,11)];
% Scale the velocity vector for better visualization
velocity_scale = 1000;
scaled_velocity = final_velocity * velocity_scale;
% Plot the position vector (as a point)
plot(final_position(1), final_position(2), 'ko', 'MarkerSize', 8, 'DisplayName', 'Final Position');
% Plot the velocity vector (as an arrow)
quiver(final_position(1), final_position(2), scaled_velocity(1), scaled_velocity(2), 'r', 'LineWidth', 2, 'DisplayName', 'Final Velocity');
% Plot the initial velocity vector with respect to the final position
initial_position_for_velocity = [Y(1,7), Y(1,8)];
initial_velocity = [Y(1, 10), Y(1, 11)];
scaled_initial_velocity = initial_velocity * velocity_scale;
quiver(initial_position_for_velocity(1), initial_position_for_velocity(2), 'g', 'LineWidth', 2, 'DisplayName', 'Initial Velocity');
L=1.5*L;
axis square
xlim([-L L]);
ylim([-L L]);
hold off
function dy = SystemsOfEquations(t,y)
global m G n mu_moon
% For each object there 3 position, and 3 velocity components. 6 in total.
dy = zeros(12,1);    % Column vector to initialize the .
    for i=1:1:n
        toplam=[0 0 0];
        for j=1:1:n
            if i==j, continue, end
            r=[ y((j-1)*6+1)-y((i-1)*6+1)  y((j-1)*6+2)-y((i-1)*6+2) y((j-1)*6+3)-y((i-1)*6+3)];
            toplam=toplam+(G*m(i)*m(j)*r/norm(r)^3);
        end
dy((i-1)*6+1) = y((i-1)*6+4);   % x coordinate of the object
dy((i-1)*6+2) = y((i-1)*6+5);   % y coordinate of the object
dy((i-1)*6+3) = y((i-1)*6+6);   % z coordinate of the object
dy((i-1)*6+4) = toplam(1)/m(i); % Velocity component along x
dy((i-1)*6+5) = toplam(2)/m(i); % Velocity component along y
dy((i-1)*6+6) = toplam(3)/m(i); % Velocity component along z
    end
end
function [value, isterminal, direction] = liftoff_start(t,y)
% This event triggers shortly after the start to initiate the ascent loop.
% It doesn't necessarily correspond to reaching a specific altitude.
value = t - 1; % Trigger after 1 second (adjust as needed)
isterminal = 1; % Continue the simulation after this event
direction = 1; % Trigger when value crosses zero from negative to positive
end
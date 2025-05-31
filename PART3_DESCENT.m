% LUNAR MISSION TRAJECTORIES | TERM PROJECT | K.G. ULUTAS
% PART 3 - POWERED DESCENT INITIATION
clear, clc, close all;
global m G n mu_moon Rm

%CONSTANTS
n=2;               % Number of the objects;
G = 6.6742867e-11; % Universal Gravitational Constant [N*m*m/kg/kg]
Re=6.3781366e6;  % Earth's radius
m2=5.9721424e24; %  Earth's mass
m1=7.3457576e22; %  Moon's mass
m3=1000;         %  Spacecraft's mass
mp = 11816.793;
mt = m3 + mp;
Rm = 1737400;
mu_earth=G*m2;
mu_moon=G*m1;
v_park_moon = sqrt(mu_moon/2097400); % Corrected parking orbit velocity
Isp = 336; % Monomethyl Hydrazin
g0 = 9.80665;
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
    b02=[1.668349469778299e+06 1.271176529678822e+06 0];        % Initial position for m2
    v02=[0.954170844583487*9.201844747344459e+02 -0.954170844583487*1.220964939316312e+03 0];       % Initial velocity for m2
    
% Initial Conditions Vector
b0=[b01 v01 b02 v02];
% Calculate the delta V before applying rocket motion

    r_sc = norm(b02(1:3));
    ra = r_sc;
    rp = Rm+15240;
    e_tr = (ra - rp) / (ra + rp) ;
    h_tr = sqrt(mu_moon*ra*(1+cosd(180)*e_tr));
    delta_v_transfer = v_park_moon - norm(v02(1:3));
    % You might need a delta-v for orbit insertion at the Moon as well.

% ÇÖZÜM ARALIĞI
tspan1=2*pi/w;
tspan=[0 tspan1]; % [s]
options = odeset('RelTol',1e-9,'AbsTol', 1e-9,'Events',@moon_altitude_descent);
[t1,y1,te,ye,ie] = ode45(@SystemsOfEquationsDescent,tspan,b0,options);
if ~isempty(te)
% Initialize from the moment of descent initiation
b0_current = ye; % Start from the event state
T = t1;
Y = y1;
descent_complete = false; % Flag to indicate descent completion

% Thrust Parameters
T_max = 15000; % Maximum thrust in Newtons
Isp = 336; % Specific impulse in seconds
g0 = 9.80665; % Standard gravity
mdot_max = T_max / (Isp * g0); % Maximum mass flow rate

for c = 1:1:700
    dmdt = 10100 / (c + 10122689 / 5000);
    m_t = -(10100 * log(abs(c + 10122689 / 5000)) - 86992.277);
    mt = mt - dmdt;
    % Compute position and velocity vectors
    r_sc_vector = b0_current(7:9);
    vel_sc = b0_current(10:12);
    fprintf('The spacecrafts initial velocity is: %d\n',norm(vel_sc));
    vel_mag = norm(vel_sc);
    vel_u = vel_sc / vel_mag;
    current_altitude = norm(r_sc_vector) - Rm;

    h = cross(r_sc_vector, vel_sc);

    gamma = acosd( dot(-r_sc_vector, vel_sc) / (norm(r_sc_vector)*norm(vel_sc))) ; % The angle between the velocity and the position vector of the spacecraft
    phi = 90 - gamma;
    
    % With respect to the moon centered inertial

    nu = -r_sc_vector / norm(r_sc_vector); %The lunar nadir vector
    vr = vel_mag * cosd(gamma) .* nu; %Radial velocity vector of the spacecraft
    vr_u = vr / norm(vr); %Unit vector of the radial velocity vector
    v_t = vel_sc - vr; %Tangential velocity vector of the spacecraft
    v_t_u = v_t / norm(v_t); %Unit vector of the tangential velocity vector of the spacecraft


    % Thrust Vector Calculation (Tangential and Radial Components)
    % This is a simplified thrust profile. You might need a more sophisticated control system.
    thrust_tangential = v_t_u * T_max; % Thrust opposite to tangential velocity
    thrust_radial = nu * T_max; % Small radial thrust for altitude control
    
    % Calculate thrust acceleration
    thrust_accel_tan = 4*thrust_tangential / mt; % Acceleration due to tangential thrust
    thrust_accel_rad = 1.22674*thrust_radial / mt; % Acceleration due to radial thrust

    % Direction of braking (opposite to velocity vector)
    if norm(vel_sc) > 1e-6
        % Apply delta-v to velocity
        b0_current(10:12) = b0_current(10:12) - v_t - thrust_accel_rad;
        fprintf("The delta v applied: v+= %f, vr= %f \n",norm(v_t), norm(thrust_accel_rad));
    else
        fprintf('Velocity is near zero, no delta-v applied.\n');
    end

    % Integrate for 1 second
    options_step = odeset('RelTol',1e-9,'AbsTol', 1e-9);
    [t_step, y_step] = ode45(@SystemsOfEquationsDescent, [0 1], b0_current, options_step);

    % Update current state
    b0_current = y_step(end, :)';
    
    % Append results (shift time by previous max)
    T = [T; T(end) + t_step(2:end)];
    Y = [Y; y_step(2:end, :)];

     % Check for landing condition
    r_sc_vector_new = b0_current(7:9);
    current_altitude_new = norm(r_sc_vector_new) - Rm;
        if norm(b0_current(10:12)) > 3 && current_altitude_new <= 0
        fprintf('Spacecraft crashed on the Moon.\n');
        fprintf('The landing velocity is: %d\n',norm(b0_current(10:12)));
        descent_complete = true;
        break; % Stop the descent loop
        else
        end
        if norm(b0_current(10:12)) <= 3 && current_altitude_new <= 0
        fprintf('Spacecraft landed on the Moon.\n');
        fprintf('The landing velocity is: %d\n',norm(b0_current(10:12)));
        descent_complete = true;
        sc_final_pos = b0_current(7:9)';
        break; % Stop the descent loop
        end
end

if ~descent_complete
    fprintf('Spacecraft did not land within %d seconds.\n',c);
    fprintf('Current velocity is: %d\n',norm(b0_current(10:12)));
end
end

position = [Y(length(T),7) Y(length(T),8)];
velocity = [Y(length(T),10) Y(length(T),11)];
Vector = [position, velocity];

%DO NOT CHANGE ANY CODE FROM THE PLOT SECTION

figure
circle = 0:pi/5000:2*pi;
hold on
plot(Y(:,1),Y(:,2),'r',...
Y(:,7),Y(:,8),'b',...
Y(length(T),1)+1737400*cos(circle),Y(length(T),2)+1737400*sin(circle),'r',...
Y(length(T),1)+2097400*cos(circle),Y(length(T),2)+2097400*sin(circle),'r--',...
Y(length(T),1)+1752640*cos(circle),Y(length(T),2)+1752640*sin(circle),'g--');

% Get the final position and velocity vectors.
final_position = [Y(end,7), Y(end,8)];
final_velocity = [Y(end,10), Y(end,11)];

% Scale the velocity vector for better visualization
velocity_scale = 1000;  % Adjust this value as needed to make the arrow visible
scaled_velocity = final_velocity * velocity_scale;

% Plot the position vector (as a point)
plot(final_position(1), final_position(2), 'ko', 'MarkerSize', 8, 'DisplayName', 'Final Position');

% Plot the velocity vector (as an arrow)
quiver(final_position(1), final_position(2), scaled_velocity(1), scaled_velocity(2), 'r', 'LineWidth', 2, 'DisplayName', 'Scaled Velocity');

% Plot the initial velocity vector with respect to the final position
initial_position_for_velocity = final_position;  % Start at the final position
initial_velocity = [Y(1, 10), Y(1, 11)];
scaled_initial_velocity = initial_velocity * velocity_scale;
quiver(initial_position_for_velocity(1), initial_position_for_velocity(2), -scaled_initial_velocity(1), -scaled_initial_velocity(2), 'g', 'LineWidth', 2, 'DisplayName', 'Initial Velocity');

L=1.5*L;
axis square
xlim([-L L]);
ylim([-L L]);

hold off

function [value, isterminal, direction] = moon_altitude_descent(t,y)
    mu_moon_local = 4902769.2;
    Rm_local = 1737400;
    r_sc = sqrt( y(7)^2 + y(8)^2 + y(9)^2 );
    ra = r_sc;
    rp_local = Rm_local + 15240;
    e_tr = (ra - rp_local) / (ra + rp_local) ;
    h_tr = sqrt(mu_moon_local*ra*(1+cosd(180)*e_tr));
    tol = 20;
    target_altitude = rp_local;
    stop_altitude_local = Rm_local;
    value_descent = abs(target_altitude - r_sc) - 100;
    isterminal_descent = 1;
    direction_descent = -1;



    value_stop = abs(stop_altitude_local - r_sc) ;
    isterminal_stop = 1;
    direction_stop = 1;
    value = value_descent; % Trigger only descent initiation
    isterminal = isterminal_descent;
    direction = direction_descent;
end

function dy = SystemsOfEquationsDescent(x,y)
global m G n

% For each object there 3 position, and 3 velocity components. 6 in total.

dy = zeros(6,1);    % Column vector to initialize the .

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

%NOT CHANGED


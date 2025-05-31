% LUNAR MISSION TRAJECTORIES | TERM PROJECT | K.G. ULUTAS

% PART 1 - 8 FIGURE TRAJECTORY 

clear, clc, close all;
global m G n

%CONSTANTS
n=3;                % Number of the objects; 
G = 6.6742867e-11;  % Universal Gravitational Constant [N*m*m/kg/kg]
Re=6.3781366e6;     % Earth's radius
Rm=1.7374e6;        %Moon radius
m1=5.9721424e24;    %  Earth's mass
m2=7.3457576e22;    %  Moon's mass
m3=1000;            %  Spacecraft's dry mass
mp = 52768.013;          % Propellant Mass
mt = m3 + mp;       % Spacecraft's total mass
mu_earth=G*m1;
mu_moon=G*m2;
v_park_earth=sqrt(mu_earth/(Re+560000));
v_park_moon=sqrt(mu_moon/(2097400));
g0 = 9.80665;

Isp = 421; % TLI Propellant J-2 Liquid

T_sc_earth = 2*pi / sqrt(mu_earth / (Re+560000)^3);
T_sc_moon = 2*pi / sqrt(mu_moon / (Re+360000)^3);

%m1=2e20; %kg
%m2=2e20; %kg
%m3=1e20; %kg

L=3.844e8; %Earth to Moon distance in m

R_soi = L*((m2/m1)^(2/5));

m=[m1 m2 m3]; % Mass vector

Lm1=m2/(m1+m2)*L;
Lm2=m1/(m1+m2)*L;
w=sqrt(G/L^3*(m1+m2+mt*(1+m1/m2)^3));
T_planet = 2*pi / w;
   
    b01=[-Lm1 0 0];        % Initial position for m1 (Earth)
%          1  2 3

    v01=[0 -w*Lm1 0];      % Initial velocity for m1 
%        4   5    6

    b02=[Lm2 0 0];        % Initial position for m2  (Moon)
%        7   8 9

    v02=[0 w*Lm2 0];      % Initial velocity for m2 
%       10  11   12
        
    %Lm2*1e-7
    b03=[-Lm1 -(Re+560000) 0];       % Initial position for m3     (Spacecraft)
%         13     14        15

    v03=[10650*cosd(20.734172) -10650*sind(20.734172) 0];       % Initial velocity for m3 
%           16                          17            18

    
% Initial Conditions Vector

b0=[b01 v01 b02 v02 b03 v03];

r_earth_sc0 = b0(13:15) - b0(1:3);

e_mech = (10650^2) / 2 - mu_earth / 6938136.6;
a = -mu_earth / e_mech;


% ÇÖZÜM ARALIĞI
tspan1=T_planet; %The moment when the spacecraft enter the Moon's circular orbit
tspan=[0 tspan1]; % [s] 
options = odeset('RelTol',1e-9,'AbsTol',1e-9,'Events',@moon_altitude_event);
[t1,y1,te,ye,ie] = ode45(@SystemsOfEquations,[0 T_planet],b0,options);
%(39.44/360)-143

if ~isempty(te)
    % Calculate required velocity for circular orbit
    r_moon_sc = [ye(13)-ye(7), ye(14)-ye(8), ye(15)-ye(9)];
    distance = norm(r_moon_sc);
    v_moon = [ye(10) ye(11) ye(12)];
    v_sc = [ye(16) ye(17) ye(18)];
    
    % Velocity relative to Moon
    v_rel_moon = v_sc - v_moon;
    v_rel_mag = norm(v_rel_moon);
    
    % Required circular velocity
    v_circular = sqrt(mu_moon/distance);
    
    % Scale velocity vector to match circular orbit speed
    v_rel_new = v_rel_moon * (v_circular/v_rel_mag);
    v_sc_new = v_moon + v_rel_new;
    
    % Update initial conditions for second integration
    b0_new = ye;
    b0_new(16:18) = v_sc_new;
    
    % Second integration - after circularization
    [t2,y2] = ode45(@SystemsOfEquations,[te T_planet],b0_new,options);
    norm(v_sc)
    % Combine results
    T = [t1; t2];
    Y = [y1; y2];
else
    % If we didn't reach the target altitude, just use the first integration
    T = t1;
    Y = y1;
    
end


figure 
circle = 0:pi/50:2*pi;

hold on
plot(Y(:,1),Y(:,2),'r'...
    ,Y(:,7),Y(:,8),'b'...
    ,Y(:,13),Y(:,14),'g'...
    , Y(length(T),7)+1737400*cos(circle),Y(length(T),8)+1737400*sin(circle),'r'...
    ,Y(length(T),7)+2097400*cos(circle),Y(length(T),8)+2097400*sin(circle),'r--'...
    ,Y(length(T),13)+50000*cos(circle),Y(length(T),14)+50000*sin(circle),'k'...
    ,Y(length(T),1)+6371136.6*cos(circle),Y(length(T),2)+6371136.6*sin(circle),'b');

L=2*L;
axis square
xlim([-L L]);
ylim([-L L]);


for i = 1:length(T)

    %Uzay aracının Dünya'ya göre yüksekliği%

    r_dunya(i,1) = sqrt( (Y(i,1) - Y(i,13))^2 + (Y(i,2) - Y(i,14))^2) - Re ;

    %Uzay aracının Ay'a göre yüksekliği%

    r_ay(i,1) = sqrt( (Y(i,7) - Y(i,13))^2 + (Y(i,8) - Y(i,14))^2 ) - 1737400 ;

    %Uzay aracının origine göre açısı

    angle_spacecraft(i,1) = atand( Y(i,14)/Y(i,13) );

    %Ay'ın origine göre açısı

    angle_moon(i,1) = atand( Y(i,8)/Y(i,7) );

    %Uzay aracının hızı

    sc_velocity_V(i,1) = sqrt( Y(i,16)^2 + Y(i,17)); 
end

r_ay(length(T),1);

sc_velocity = sqrt( Y(length(T),16)^2 + Y(length(T),17)^2);

% fprintf("Spacecraft Velocity Vector %d i + %d j",Y(length(T),16),Y(length(T),17));

function [value, isterminal, direction] = moon_altitude_event(t, y)
    Rm = 1.7374e6;
    target_altitude = 360000; % meters
    tolerance = 46;           % meters

    % Vector from Moon to spacecraft
    r_moon_sc = y(13:15) - y(7:9);
    current_altitude = norm(r_moon_sc) - Rm;

    value = abs(current_altitude - target_altitude) - tolerance;
    isterminal = 1;
    direction = -1;

    % G = 6.6742867e-11;  % Universal Gravitational Constant [N*m*m/kg/kg]
    % Mm=7.3457576e22;
    % T_co = 2*pi*sqrt( (Rm+360e3)^3 / (G*Mm)); % Period of the spacecraft around the moon in circular orbit
    % tolerance = 29;         % Seconds
    % value_co = abs(t - 2.581131994659712e+05 + T_co) - tolerance; % How to calculate the time passed until the spacecraft reahes the circular orbit and make delta V?
    % isterminal_co = 1;
    % direction_co = 0;
    % value_co
    % value = [value_al, value_co];
    % isterminal = [isterminal_al, isterminal_co];
    % direction = [direction_al, direction_co];
    % 

end

deltav = norm(v03 - [v_park_earth 0 0]);
m_new = mt - mt*(1-exp(1)^(-deltav/(Isp*g0)));

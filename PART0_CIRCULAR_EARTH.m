clear, clc, close all;
global m G n

%CONSTANTS
n=3;               % Number of the objects; 
G = 6.6742867e-11; % Universal Gravitational Constant [N*m*m/kg/kg]
Re=6.3781366e6; % Earth's radius
Rm=1.7374e6;   %Moon radius
m1=5.9721424e24; %  Earth's mass
m2=7.3457576e22; %  Moon's mass
m3=1000;         %  Spacecraft's mass
mu_earth=G*m1;
mu_moon=G*m2;
v_park_earth=sqrt(mu_earth/(Re+560000));
v_park_moon=sqrt(mu_moon/(2097400));
N = 6000;

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
w=sqrt(G/L^3*(m1+m2+m3*(1+m1/m2)^3));
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

    v03=[7.5795966e3 0 0];       % Initial velocity for m3 
%           16                          17            18
    
% Initial Conditions Vector

b0=[b01 v01 b02 v02 b03 v03];

% ÇÖZÜM ARALIÐI
tspan1=2*pi*sqrt(6938136.6^3/mu_earth);
tspan=[0 tspan1]; % [s] 
options = odeset('RelTol',1e-9,'AbsTol', 1e-9); 
[T,Y] = ode45(@SystemsOfEquations,tspan,b0,options);

circle = 0:0.001*pi:2*pi;
figure 

hold on
plot(Y(:,1),Y(:,2),'r',Y(:,7),Y(:,8),'b',Y(:,13),Y(:,14),'g'...
    ,Y(length(T),13)+50000*cos(circle),Y(length(T),14)+50000*sin(circle),'k'...
    ,Y(length(T),1)+6371136.6*cos(circle),Y(length(T),2)+6371136.6*sin(circle),'b');

L=2*L;
axis square
xlim([-L L]);
ylim([-L L]);


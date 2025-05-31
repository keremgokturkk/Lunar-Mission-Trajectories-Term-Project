% LUNAR MISSION TRAJECTORIES | TERM PROJECT | K.G. ULUTAS

% PART 2 - CIRCULAR ORBIT AROUND THE MOON

clear, clc, close all;
global m G n

%CONSTANTS
n=2;               % Number of thsqrt( Y(i,10)^2 + Y(i,11)^2)e objects; 
G = 6.6742867e-11; % Universal Gravitational Constant [N*m*m/kg/kg]
Re=6.3781366e6;  % Earth's radius
m2=5.9721424e24; %  Earth's mass
m1=7.3457576e22; %  Moon's mass
m3=1000;         %  Spacecraft's mass
mp = 17258.372;
mt = m3 + mp;

mu_earth=G*m2;
mu_moon=G*m1;

% t=0:0.01:2*pi;     % 

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
    v02=[9.201844747344459e+02 -1.220964939316312e+03 0];       % Initial velocity for m2 
    
% Initial Conditions Vector

b0=[b01 v01 b02 v02];

% ÇÖZÜM ARALIÐI
tspan1=2*pi/w;
tspan=[0 tspan1]; % [s] 
options = odeset('RelTol',1e-9,'AbsTol', 1e-9); 
[T,Y] = ode45(@SystemsOfEquations,tspan,b0,options);


figure 
circle = 0:pi/50:2*pi;
hold on
plot(Y(:,1),Y(:,2),'r',Y(:,7),Y(:,8),'b',Y(length(T),1)+1737400*cos(circle),Y(length(T),2)+1737400*sin(circle),'r',Y(length(T),1)+2097400*cos(circle),Y(length(T),2)+2097400*sin(circle),'r--');

L=1.5*L;
axis square
xlim([-L L]);
ylim([-L L]);

Isp = 336;
g0 = 9.80665;
v_park_moon = sqrt(4902769.2 / 2097400);
deltav = v_park_moon - norm(v02);
m_new = mt - mt*(1-exp(1)^(deltav/(Isp*g0)));




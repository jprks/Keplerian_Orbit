% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
    Title: Keplerian trajectory determination for Re-Entry Vehicles
    Author: James Emerson Parkus
    Date: 07/10/2018
    Purpose: This script will determine the orbit for a re-entry vehicles
    given a satellite mass, initial altitude, entry angle, and arguement of
    apoapsis. 
    
    Output: entry velocity, Delta V to de-orbit from LEO

    G            - gravitational constant [m^3/kg*s^2]
    R_e          - earth radius [m]
    M            - earth (planet) mass [kg]
    mu           - gravitational parameter [m^3/s^2]
    m            - satellite mass [kg]
    alt          - altitude [m]
    entry angle  - entry angle of RV to atmosphere [deg]
    r_a          - arguement of apoapsis [m]
    r_atmos      - radius of atmosphere [m]
    a            - semimajor axis [m]
    e            - eccentricity [-]
    theta_i      - intersect angle [deg]
    gamma        - intersect flight angle [deg]
    orbit_geo    - geometry of satellite orbit [m]
    h            - angular momentum [m^2/s]
    flight_path  - flight path of entire desired flight [deg]
    v_r          - radial velocity [m/s]
    v_p          - perpendicular velocity [m/s]
    v_orbit      - orbital velocity [m/s]   
    y            - State vector for RV [m m/s]
%}
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clc
close all
clear

%% Constants
G = 6.6742*10^(-11); % m^3/kg*s^2
R_e = 6.378*10^6; % m - Earth Radius
M = 5.97*10^24; % kg - Earth mass
mu = G*M; % m^3/s^2 - Gravitational parameter

%% Inputs
m = 4; % kg - Satellite mass
alt = 400000; % m - Altitude of r0
entry_angle = deg2rad(5); % deg - Entry Angle
r_a = R_e + alt; % m - Arguement of Apoapsis
r_atmos = R_e + 100000; % m - Radius of Atmosphere

%% Allocating Memory
n = 2*pi/deg2rad(0.1);
p = R_e/1000;

r = zeros(1,n);
earth_theta = zeros(1,n);
theta_int = zeros(1,p);
gamma = zeros(1,p);
r_periapsis = zeros(1,p);
flight_path = zeros(1,n);
v_r = zeros(1,n);
v_p = zeros(1,n);
alpha = zeros(1,n);
atmos = zeros(1,n);
r_ellipse = zeros(1,n);
a = zeros(1,p);
v_orbit_mag = zeros(1,n);
e = zeros(1,n);

%% Earth
l = 1;
for theta = 0:deg2rad(0.1):2*pi
    r(l) = R_e;
    atmos(l) = r_atmos;
    earth_theta(l) = theta;
    l = l + 1;
end

%% Intercept Conditions
i = 1; % Iteration counter
for r_p = 0:1000:R_e
    a(i) = 1/2*(r_p+r_a); % m - Semimajor Axis
    e(i) = (r_a - r_p)/(r_a + r_p); % Eccentricity
    theta_int(i) = acos(1/e(i)*(a(i)*(1 - e(i)^2)/r_atmos - 1)); % Intersect angle
    gamma(i) = atan(e(i)*sin(theta_int(i))/(1 + e(i)*cos(theta_int(i)))); % Intersect flight path angle
    r_periapsis(i) = r_p;
    i = i + 1;
end

for t = 2:1:p
    if gamma(t-1) > entry_angle && gamma(t+1) < entry_angle
        entry_flight_path = gamma(t);
        entry_intersect_angle = theta_int(t);
        u = t;
    end
end

%% Analysis
r_periapsis = r_periapsis(u);
a = a(u);
e = e(u);

h = sqrt(mu*a*(1 - e^2));

i = 1;
for angle = 0:deg2rad(0.1):2*pi
    r_ellipse(i) = a*(1 - e^2)/(1 + e*cos(angle));
    flight_path(i) = atan(e*sin(angle)/(1 + e*cos(angle)));
    v_r(i) = mu/h*e*sin(angle);
    v_p(i) = mu/h*(1 + e*cos(angle));
    alpha(i) = angle;
    i = i + 1;
end

v_r = v_r';
v_p = v_p';

v_orbit = [v_r v_p];

for j = 1:1:n+1
    v_orbit_mag(j) = sqrt(v_orbit(j,1)^2+v_orbit(j,2)^2);
end

alpha = alpha';
r_ellipse = r_ellipse';
r = r';
v_orbit_mag = v_orbit_mag';
atmos = atmos';

entry_velocity = sqrt((mu/h*e*sin(entry_intersect_angle))^2+(mu/h*(1 + e*cos(entry_intersect_angle)))^2);

y = [r_ellipse v_orbit_mag];

v_apoapsis = min(v_orbit_mag);

v_circle = sqrt(mu/r_a);

delta_v = v_circle - v_apoapsis;
%% Conversions
entry_angle = rad2deg(entry_angle);
flight_path = rad2deg(flight_path);
theta_int = rad2deg(theta_int);
gamma = rad2deg(gamma);
entry_flight_path = rad2deg(entry_flight_path);

%% Output

fprintf('\nThe entry velocity is %g m/s with an angle of %g deg from the local horizon.\n',entry_velocity,entry_angle)
fprintf('\nThe necessary change in velocity to de-orbit the satellite is %g m/s.\n',delta_v)

hold on

polar(alpha,r,'g');
polar(alpha,atmos,'b');
polar(alpha,r_ellipse,'k');
title('RV Orbit');
xlabel('m');
ylabel('m');
legend('Earth','Atmosphere','RV Orbit')

axis equal

hold off

figure()
plot(alpha,v_orbit_mag);
title('Orbital Velocity over flight');

figure()
plot(alpha,flight_path);
title('Flight path angle over flight');
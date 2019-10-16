% Given the wheel speeds (in rpm) and launch angle, this
% program computes the ball's trajectory, accounting for the effects
% of gravity, air drag, and the Magnus effect.

% Wheel names:
% x1 and x2 are the wheels along the x axis, where x1 is on the negative x
% side of the ball and x2 is on the positive x side.
% Similar definition for wheels z1 and z2.

% Coordinates: 
% x and y are in the horizontal plane, z is height above ground
% The machine is located at (x=0, y=0).
% y points in the same direction as the machine, x points "to the right"

%%% INITIAL VALUES (SET BY USER) %%%

% Typical curving kick: velocity = 21.1 m/s, pitch = 26.2, launch spin rate of 100 rad/s

% Launch position
x = 0; % meters "to the right" of the ball launcher
y = 0; % meters "straight ahead" of the ball launcher
z = 1; % height (meters) above ground

% Speed, rotation and launch angle
rpm_x1 = 8e3;
rpm_x2 = 6e3;
rpm_z1 = 7e3;
rpm_z2 = 7e3;

yaw = 0; % degrees left or right (right is positive)
pitch = 25; % degrees; angle above the horizontal



%%% SIMULATION PARAMETERS %%%

mps2mph = 2.23694;

% Properties of tennis ball
d = 0.0670; % Diameter (meters); International Tennis Federation (ITF) defines the official diameter as 6.54–6.86 cm (2.57–2.70 inches)
m = 0.0577; % Mass of ball (kg); Balls must have masses in the range 56.0–59.4 g (1.98–2.10 ounces).

% Properties of air (at 20 deg C)
rho = 1.20; % Density (kg/m^3)
mu = 1.82e-5; % Dynamic viscosity (N*s/m^2)

% Wheel dimensions
d_wheel = 0.0698;
r_wheel = d_wheel / 2;

% Calculate each wheel's required rpm 
r_ball = d / 2;

wheel_rpm = [rpm_x1 rpm_x2 rpm_z1 rpm_z2];
wheel_w = wheel_rpm * 2*pi / 60; % rad/s
wheel_v = wheel_w*r_wheel; % velocity at edge of wheel [m/s]

V = mean([wheel_v(1) wheel_v(2)]);
V_mph = V * mps2mph;
% if V ~= mean([wheel_v(3) wheel_v(4)])
%     error('Average speed in x-direction must equal average speed in z-direction!')
% else V_mph = V * mps2mph;
% end
topspin = (wheel_v(4) - V) / r_ball;
sidespin = (wheel_v(2) - V) / r_ball;

% Initial rotation (assume this will stay unchanged during the ball's flight [suggested by papers on soccer ball])
omega = [-topspin, 0, sidespin]; % initial angular velocity [rad/s]



%%% START SIMULATION %%%

% Simulation parameters
ts = 0.001; % timestep length [sec]

% Initial calculations: Vectorize the inputs
X = [x y z];

Vx = V*cosd(pitch)*sind(yaw);
Vy = V*cosd(pitch)*cosd(yaw);
Vz = V*sind(pitch);
v = [Vx Vy Vz]; % initial velocity [m/s]

Omega = norm(omega);



A = pi*(d/2)^2; % Cross-sectional area of ball

t = 0;
a = [0 0 0];
while X(end,3) > 0 % height above ground > 0
    my_v = v(end,:);
    V = norm(my_v);
	n_v = my_v / V; % unit vector for velocity
	
    Re = rho*V*d/mu; % Reynold's number
	SP = Omega * (d/2) / (V); % "Spin Parameter": ratio of (maximum circumferential velocity on ball surface) / (free-stream velocity)
	
    Cd = 0.65; % may range from 0.625 to 0.675
    Cm = 1/(2.022+0.981 / SP); % +/- 0.05
    
	F = 0.5*rho*A*V^2;
	F_d = F*Cd*-n_v; % Drag force vector; direction is opposite to velocity
	
    if Omega > 0
        n_omega = omega / Omega;
    else n_omega = [0 0 0];
    end
	n_m = cross(n_omega, n_v);
	F_m = F*Cm*n_m; % Magnus force vector; direction is 
    
	F_g = [0 0 -9.8*m];
	
	F_tot = F_g + F_d + F_m;
    
    % get a, v, and X
    a_res = F_tot / m; % resultant acceleration [m/s^2]
	a = [a; a_res];
    v_res = my_v + a_res*ts;
    v = [v; v_res]; % Get velocity vector for next timestep by integrating acceleration
    X_res = X(end,:) + v_res*ts;
    X = [X; X_res]; % Get position for next timestep by integrating velocity
    
    t = t + ts;
end
figure
scatter3(X(:,1), X(:,2), X(:,3))
title(['t = ' num2str(length(X) / 1000) ' seconds'])
xlabel('x (meters)')
ylabel('y (meters)')
zlabel('Height (meters)')
axis equal
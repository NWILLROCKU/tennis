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
X0 = [x y z];

yaw = 0; % degrees left or right (right is positive)


%%% SIMULATION PARAMETERS %%%

mps2mph = 2.23694;

% Properties of tennis ball
d = 0.0670; % Diameter (meters); International Tennis Federation (ITF) defines the official diameter as 6.54–6.86 cm (2.57–2.70 inches)
mass = 0.0577; % Mass of ball (kg); Balls must have masses in the range 56.0–59.4 g (1.98–2.10 ounces).
r_ball = d / 2;

% Properties of air (at 20 deg C)
rho = 1.20; % Density (kg/m^3)
mu = 1.82e-5; % Dynamic viscosity (N*s/m^2)

% Wheel dimensions
d_wheel = 0.0698; % [m]
r_wheel = d_wheel / 2;



%%% START SIMULATION %%%

% Simulation parameters
ts = 0.001; % timestep length [sec]

A = pi*(d/2)^2; % Cross-sectional area of ball

rpms_x1 = 1000:1000:10e3;
rpms_x2 = 1000:1000:10e3;
pitches = 0:5:45; % pitches to test (pitch = angle above horizontal [deg])

m = length(rpms_x1);
n = length(rpms_x2);
p = length(pitches);
dest_x=[];
dest_y=[];
% dest_x=zeros(m,n,p);
% dest_y=zeros(m,n,p);
for i = 1:m
    rpm_x1 = rpms_x1(i);
    for j = 1:n
        rpm_x2 = rpms_x2(j);
        wheel_rpm = [rpm_x1 rpm_x2];
        wheel_w = wheel_rpm * 2*pi / 60; % rad/s
        wheel_v = wheel_w*r_wheel; % velocity at edge of wheel [m/s]

        V0 = mean(wheel_v); % speed at center of ball
        wheel_v = [wheel_v V0 V0];
        topspin = (wheel_v(4) - V0) / r_ball;
        sidespin = (wheel_v(2) - V0) / r_ball;

        % Initial rotation (assume this will stay unchanged during the ball's flight [suggested by papers on soccer ball])
        omega = [-topspin, 0, sidespin]; % initial angular velocity [rad/s]
        Omega = norm(omega);
        
        for k = 1:p
            pitch = pitches(k);
            
            % Initialize X, v, a, t
            X = [x y z];

            Vx = V0*cosd(pitch)*sind(yaw);
            Vy = V0*cosd(pitch)*cosd(yaw);
            Vz = V0*sind(pitch);
            v = [Vx Vy Vz]; % initial velocity [m/s]

            a = [0 0 0];
            t = 0;
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

                F_g = [0 0 -9.8*mass];

                F_tot = F_g + F_d + F_m;

                % get a, v, and X
                a_res = F_tot / mass; % resultant acceleration [m/s^2]
                a = [a; a_res];
                v_res = my_v + a_res*ts;
                v = [v; v_res]; % Get velocity vector for next timestep by integrating acceleration
                X_res = X(end,:) + v_res*ts;
                X = [X; X_res]; % Get position for next timestep by integrating velocity

                t = t + ts;
            end
            dest_x=[dest_x; X(end, 1)];
            dest_y=[dest_y; X(end, 2)];
%             dest_x(i,j,k) = X(end, 1);
%             dest_y(i,j,k) = X(end, 2);
        end
    end
end
figure
scatter(dest_x, dest_y)
xlabel('x [m]')
ylabel('y [m]')
axis equal
grid on
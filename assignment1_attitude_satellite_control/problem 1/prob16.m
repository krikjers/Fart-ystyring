    clear all;
close all;
clc;

% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and Håkon H. Helgesen

%% USER INPUTS
h = 0.1;                     % sample time (s)
N  = 3200;                    % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;

I = diag([m*r^2 m*r^2 m*r^2]);       % inertia matrix
I_inv = inv(I);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

%%%
phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = 20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,14+3+3);        % memory allocation

kp = 20;
kd = 400;

%% FOR-END LOOP
for i = 1:N+1,
   t = (i-1)*h;                  % time
   
   
   phi_d = 0;
   theta_d = deg2rad*15*cos(0.1*t);
   psi_d = deg2rad*10*sin(0.05*t);

   q_d = euler2q(phi_d, theta_d, psi_d);
   
   q_d_bar = quatconj(q_d')';          %find quaternion conjugate - column vector is the result
   
   q_tilde = quatmultiply(q_d_bar', q')';
   
   e_tilde = q_tilde(2:4);
   
   %%%%
   
   phi_d_dot = 0;
   theta_d_dot = -15*sin(0.1*t)*0.1*deg2rad;
   psi_d_dot = 10*cos(0.05*t)*0.05*deg2rad;
   
   THETA_d = [phi_d; theta_d; psi_d];
   THETA_dot = [phi_d_dot; theta_d_dot; psi_d_dot];
   
   T_inverse = [1 0  (-sin(theta_d)); 
                0 cos(phi_d) cos(theta_d)*sin(phi_d);
                0  (-sin(phi_d)) cos(theta_d)*cos(phi_d)];
            
            
   w_d = T_inverse*THETA_dot;
   
   w_i = w;
   w_d_tilde = w_i - w_d;
   e_i_tilde = q_tilde(2:4);
   %%%
   %tau = [0.5 1 -1]';               % control law  %column vector
   tau = (-kd*w_d_tilde -kp*e_i_tilde);
   
 
   
   

   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau' w_d' THETA_d'];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
   
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);

w_d     = rad2deg*table(:, 15:17);
THETA_d   = rad2deg*table(:, 18: 20); 


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');






figure (4); clf;
hold on;
plot(t, w_d(:,1), 'b');
plot(t, w_d(:,2), 'r');
plot(t, w_d(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('w_d');
xlabel('time [s]'); 
ylabel('deg');

figure (5); clf;
hold on;
plot(t, THETA_d(:,1), 'b--');
plot(t, THETA_d(:,2), 'r--');
plot(t, THETA_d(:,3), 'g--');
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
hold off;
grid on;
legend('\phi_d', '\theta_d', '\psi_d', '\phi', '\theta', '\psi');
title('Desired Euler angles vs actual Euler angles');
xlabel('time [s]'); 
ylabel('deg');


figure (6); clf;
hold on;
plot(t, w_d(:,1), 'b--');
plot(t, w_d(:,2), 'r--');
plot(t, w_d(:,3), 'g--');
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('x_d', 'y_d', 'z_d', 'x', 'y', 'z');
title('\omega desired vs actual \omega');
xlabel('time [s]'); 
ylabel('deg');
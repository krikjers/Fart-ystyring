clear all;

%% 1.5
t = 3;
   phi = 0;
   theta = 15*cos(0.1*t);
   psi = 10*sin(0.05*t);
   q_d = euler2q(phi, theta, psi);
   
   q_d_bar = quatconj(q_d')';          %find quaternion conjugate - column vector is the result
   q = [1;2;3;4];
   q_tilde = quatmultiply(q_d_bar', q')';
   
   e_tilde = q_tilde(2:4);
   
   
   %% 1.6
   clear all;
   t = 80;
   phi_d = 0;
   theta_d = 15*cos(0.1*t);
   psi_d = 10*sin(0.05*t);
   
   phi_d_dot = 0;
   theta_d_dot = -15*sin(0.1*t)*0.1;
   psi_d_dot = 10*cos(0.05*t)*0.05;
   
   THETA = [phi_d; theta_d; psi_d];
   THETA_dot = [phi_d_dot; theta_d_dot; psi_d_dot];
   
   T_inverse = [1 0  (-sin(theta_d)); 
                0 cos(phi_d) cos(theta_d)*sin(phi_d);
                0  (-sin(phi_d)) cos(theta_d)*cos(phi_d)];
            
            
   w_d = T_inverse*THETA_dot
            
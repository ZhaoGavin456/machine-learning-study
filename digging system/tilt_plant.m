function state = tilt_plant(u, F_T, x, dt)
% parameters
n_tilt = 1; % number of tilt cylinder;
m = 30*n_tilt;
k_pist = 5000; % piston spring coefficient  
b_pist = 28000; % piston damping coefficient  
tilt_A = 1 * (pi*(165.1e-3)^2/4); % tilt piston A area, m^2;
tilt_B = 1 * (pi*(165.1e-3)^2/4 - pi*(88.9e-3)^2/4); % tilt piston B area, m^2;
tilt_ave = 0.5*(tilt_A+tilt_B); % average area
tilt_delta = 0.5*(tilt_A-tilt_B);  % area difference
A_p = tilt_ave; % pressure area of the piston
A_spool = pi*5e-3*1e-3; % spool area
storke_tilt = 0.6; % piston stroke
beta = 1.8e9; % bulk modulus
V_t = storke_tilt*tilt_ave; % volume for tilt chamber
C_d = 0.6; % discharge coefficient
rho = 840; % density
P_s = 35e6; % supply pump pressure

% combination of parameters; 
% reference paper: High-Gain-Observer-Based Integral Sliding Mode Control 
% for Position Tracking of Electrohydraulic Servo Systems
% equation 6, 7; 
% I ignore leakage, C_tl = 0; Hence, my h_2 is h_3 in the paper. 
a1 = k_pist/m;
a2 = b_pist/m;
a3 = A_p/m;
h1 = beta*(tilt_A^2+tilt_B^2)/V_t/tilt_ave;
h2 = beta*(tilt_A^2+tilt_B^2)/V_t/tilt_ave*C_d*A_spool*sqrt(2/rho/(tilt_A^3+tilt_B^3));

A = [0, 1, 0;
    0, 0, a3;
    0, 0, 0];
A1x = [0;
    -a1*x(1)-a2*x(2)+a3*x(3)-F_T/m;
    -h1*x(2) ;];

phixu = [0;0; h2*sqrt(real(tilt_ave*P_s-sign(u)*(tilt_delta*P_s-tilt_ave*x(3))))*u;];

x_dot = A*x + A1x + phixu; % state dot

state = x_dot*dt + x; 

end
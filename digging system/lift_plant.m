function state = lift_plant(u, F_L, x, dt)
% parameters
n_lift = 2; % number of lift piston;
m = 30*n_lift;  % ----------- mass of lift piston
k_pist = 5000;  % piston spring coefficient  
b_pist = 28000/1; % piston damping coefficient
lift_A = n_lift * (pi*(133.4e-3)^2/4); % lift piston A area, m^2;
lift_B = n_lift * (pi*(133.4e-3)^2/4 - pi*(76.2e-3)^2/4); % lift piston B area, m^2;
lift_ave = 0.5*(lift_A+lift_B);  % average area
lift_delta = 0.5*(lift_A-lift_B); % area difference
A_p = lift_ave; % pressure area of the piston
A_spool = pi*5e-3*1e-3;  % ----------- spool area
beta = 1.8e9; % bulk modulus
storke_lift = 0.845; % piston stroke
V_t = storke_lift*lift_ave;  % ----------- volume for lift chamber
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
h1 = beta*(lift_A^2+lift_B^2)/V_t/lift_ave;
h2 = beta*(lift_A^2+lift_B^2)/V_t/lift_ave*C_d*A_spool*sqrt(2/rho/(lift_A^3+lift_B^3));

A = [0, 1, 0;
    0, 0, a3;
    0, 0, 0];
A1x = [0;
    -a1*x(1)-a2*x(2)+a3*x(3)-F_L/m;
    -h1*x(2) ;];

phixu = [0;0; h2*sqrt(real(lift_ave*P_s-sign(u)*(lift_delta*P_s-lift_ave*x(3))))*u;];

x_dot = A*x + A1x + phixu;  % state dot

state = x_dot*dt + x; 

end
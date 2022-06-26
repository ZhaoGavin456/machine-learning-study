clear all; close all; clc;

%% parameters
load('created_yd_4_28_2022.mat'); % references for tracking

dt = 1e-4;  % time step
t = tnew';  % transpose time array
tend = tnew(end); % end time 

LF = 6e5; % lift load force magnitude
F_lu =  LF/5*sin(t); % lift load force uncertainty 
F_lu = F_lu'; 
F_L = FL + F_lu; % lift load force 
F_L = F_L'; 

TF = 2.5e5; % tilt load force magnitude
F_tu = TF/5*sin(t); % tilt load force uncertainty
F_tu = F_tu';
F_T = FT + F_tu;  % tilt load force uncertainty 
F_T = F_T';

ydl = ydl';  % lift piston desired trajectory
ydl1 = ydl1'; % lift piston desired velocity
ydl2 = ydl2'; % lift piston desired acceleration
ydl3 = ydl3';

ydt = ydt';  % tilt piston desired trajectory
ydt1 = ydt1'; % tilt piston desired velocity
ydt2 = ydt2'; % tilt piston desired acceleration
ydt3 = ydt3';

%% system simulation

k_option_fb = 1;  % feedback (FB) is adopted (true)

% lift part ----------------------------------------------------
xl = zeros(3,tend/dt+1); % lift state
xl(:,1) = [0;0;0];  % lift initial states

u_fbl = zeros(1,tend/dt); % feedback control initial values for lift
ul = zeros(1,tend/dt); % control initial values for lift

k_pl = 10*2; % P gain
k_il = 0.1*2;  % I gain

err_l = zeros(1,tend/dt+1); % error intitial values
int_err_l = zeros(1,tend/dt+1); % integral error intitial values

% lift system update
for i = 1:1:tend/dt

    yl(1,i) = xl(1,i);
    err_l(1,i) = ydl(1,i) - yl(1,i);  % tracking error
    int_err_l(1,i+1) = err_l(1,i) + int_err_l(1,i);  % integral error
    u_fbl(1,i) = k_pl*err_l(1,i) + k_il*int_err_l(1,i);  % feedback control input
     
    ul(1,i) =  k_option_fb*u_fbl(1,i); %  FB 

    % lift control input saturation 
    if  ul(1,i) >1 
        ul(1,i) = 1;
    end
    if ul(1,i) <-1 
        ul(1,i) = -1;
    end

    xl(:,i+1) = lift_plant(ul(1,i), F_L(i), xl(:,i), dt);


end 


% tilt part ----------------------------------------------------
xt = zeros(3,tend/dt+1); % tilt state
xt(:,1) = [0;0;0];  % tilt initial states

u_fbt = zeros(1,tend/dt); % feedback control initial values for tilt
ut = zeros(1,tend/dt); % control initial values for tilt

k_pt = 10*2; % P gain
k_it = 0.1*2;  % I gain

err_t = zeros(1,tend/dt+1); % error intitial values
int_err_t = zeros(1,tend/dt+1); % integral error intitial values
    

% tilt system update
for i = 1:1:tend/dt

    yt(1,i) = xt(1,i);
    err_t(1,i) = ydt(1,i) - yt(1,i);  % tracking error
    int_err_t(1,i+1) = err_t(1,i) + int_err_t(1,i);  % integral error
    u_fbt(1,i) = k_pt*err_t(1,i) + k_it*int_err_t(1,i);  % feedback control input
     
    ut(1,i) =  k_option_fb*u_fbt(1,i); %  feedback (FB) 
    
    % lift control input saturation 
    if  ut(1,i) >1 
        ut(1,i) = 1;
    end
    if ut(1,i) <-1 
        ut(1,i) = -1;
    end
    
    xt(:,i+1) = tilt_plant(ut(1,i), F_T(i), xt(:,i), dt);


end 



%% plot results
%------------------------------------
figure()
subplot(2,1,1)
plot(t,F_L,'r',t,FL,'b',t,F_lu,'k--',"lineWidth", 1.2)
legend('Total','Nominal','Uncertainty','Location','best')
title('Lift Load Force (N)')
xlabel('time (s)')
subplot(2,1,2)
plot(t,F_T,'r',t,FT,'b',t,F_tu,'k--',"lineWidth", 1.2)
legend('Total','Nominal','Uncertainty','Location','best')
title('Tilt Load Force (N)')
xlabel('time (s)')

figure()
subplot(2,1,1) % ---- output iteration 10 vs reference
plot(t(1:end-1), yl(:)*1000, t(1:end-1), ydl(1:end-1)*1000);
title('Lift refernce v/s output (mm)');
legend('FB output','reference','location','best')
subplot(2,1,2) % ---- error
plot(t, err_l*1000);title('FB lift error [mm]');grid on
xlabel('time (s)')

figure()
subplot(2,1,1) % ---- output iteration 10 vs reference
plot(t(1:end-1), yt(:)*1000, t(1:end-1), ydt(1:end-1)*1000);
title('Tilt refernce v/s output (mm)');
legend('FB output','reference','location','best')
subplot(2,1,2) % ---- error
plot(t, err_t*1000);title('FB tilt error [mm]');grid on
xlabel('time (s)')

figure()
subplot(2,1,1)
plot(t(1:end-1), ul(1,:))
title('lift FB input');
xlabel('time (s)')
subplot(2,1,2)
plot(t(1:end-1), ut(1,:))
title('tilt FB input');
xlabel('time (s)')


%% getting bucket tip position

% combine lift and tilt trajectory to buckt tip trajectory
% addpath 'G:\My Drive\[Modeling] full order model\modeling with steering control\simulationfiles'

% these files (from TAMU collaborators) are used to convert the lift and tilt piston motion to 
% bucket motion. The vehicle motion is fixed as reference.
run('BucketDynamicsParameters.m')
run('SoilParameters.m')
load('soil_params.mat')

% bucket tip trajectory calculation with reference
[theta_t, x_t, y_t, dot_theta_t, dot_x_t, dot_y_t, ddot_theta_t, ...
    ddot_x_t, ddot_y_t] = ...
    getGeometrywithAcc(ydl+lift_ini*ones(1,length(tnew)), ...
    ydt+tilt_ini*ones(1,length(tnew)), ydl1, ...
    ydt1, ydl2, ydt2, ...
    veh_disp', y_A_init,...
    veh_vel', 0, ...
    veh_acc', 0);

% bucket tip trajectory calculation with simulated trajectory
[theta_t_a, x_t_a, y_t_a, dot_theta_t_a, dot_x_t_a, dot_y_t_a, ddot_theta_t_a, ...
    ddot_x_t_a, ddot_y_t_a] = ...
    getGeometrywithAcc(xl(1,1:end)+lift_ini*ones(1,length(tnew)), ...
    xt(1,1:end)+tilt_ini*ones(1,length(tnew)), xl(2,1:end), ...
    xt(2,1:end), ydl2, ydt2, ...
    veh_disp', y_A_init,...
    veh_vel', 0, ...
    veh_acc', 0);

%% plot bucket tip trajectory 
% bucket tip trajectory 
figure()
plot(x_t-x_t(1), y_t-y_t(1),'r', "lineWidth", 1.2)
hold on;
plot(x_t_a - x_t_a(1), y_t_a - y_t_a(1),'b--',"lineWidth", 1.2)
legend('Reference','Actual','Location','northwest')
title("FB Bucket Tip Trajectory")
xlabel({'Horizontal Distance (m)'},'FontWeight','bold')
ylabel({'Vertical Distance (m)'},'FontWeight','bold')


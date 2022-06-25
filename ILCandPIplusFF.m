clear all; close all; clc;
option=menu('Control Method to Use','1. Iterative Learning Control','2. ILC + FB','3. FF + ILC + FB','Exit');

%% parameters
load('created_yd_4_28_2022.mat');

dt = 1e-4;  % time step
n = 10; % n trails
t = tnew';
tend = tnew(end); % end time 

LF = 6e5; % lift load force magnitude
F_lu =  LF/5*sin(t); % lift load force uncertainty 
F_lu = F_lu'; 
F_L = FL + F_lu; % lift load force 
F_L = F_L'; 

TF = 2.5e5; % tilt force magnitude
F_tu = TF/5*sin(t); % tilt load force uncertainty
F_tu = F_tu';
F_T = FT + F_tu;  % tilt load force uncertainty 
F_T = F_T';

ydl = ydl';  % lift piston desired trajectory
ydl1 = ydl1'; % lift piston desired velocity
ydl2 = ydl2'; % lift piston desired acceleration
ydl3 = ydl3';

ydt = ydt';  % tilt piston desired trajectory
ydt1 = ydt1';
ydt2 = ydt2';
ydt3 = ydt3';


%% feedforward controller for lift function
beta = 1.8e9;
C_d = 0.6;
rho = 840;
P_s = 35e6;

n_lift = 2;
m_lift = 30*n_lift*0.9;
k_lift = 5000;
b_lift = 28000/1;
lift_A = n_lift * (pi*(133.4e-3)^2/4); % lift piston A area, m^2;
lift_B = n_lift * (pi*(133.4e-3)^2/4 - pi*(76.2e-3)^2/4); % lift piston B area, m^2;
lift_ave = 0.5*(lift_A+lift_B);
lift_delta = 0.5*(lift_A-lift_B); 
% A_p = lift_ave; % piston area
A_spool_l = pi*5e-3*1e-3*1.1;
storke_lift = 0.845;
V_tl = storke_lift*lift_ave*0.9;
al1 = k_lift/m_lift;
al2 = b_lift/m_lift;
al3 = lift_ave/m_lift;
hl1 = beta*(lift_A^2+lift_B^2)/V_tl/lift_ave;
hl2 = beta*(lift_A^2+lift_B^2)/V_tl/lift_ave*C_d*A_spool_l*sqrt(2/rho/(lift_A^3+lift_B^3));

u_FFdl = (1/al3*(al1*ydl1+al2*ydl2+ydl3)+hl1*ydl1);
u_FFnl = hl2*sqrt(real(lift_ave*P_s-sign(u_FFdl).*(lift_delta*P_s-lift_ave*(al1*ydl+al2*ydl1+ydl2)/al3)));
u_FFl = u_FFdl./u_FFnl;  % feedforward control input for lift


%% feedforward controller for tilt function
n_tilt = 1; % number of tilt cylinder;
m_tilt = 30*n_tilt*0.9;
k_tilt = 5000;
b_tilt = 28000;
tilt_A = 1 * (pi*(165.1e-3)^2/4); % tilt piston A area, m^2;
tilt_B = 1 * (pi*(165.1e-3)^2/4 - pi*(88.9e-3)^2/4); % tilt piston B area, m^2;
tilt_ave = 0.5*(tilt_A+tilt_B);
tilt_delta = 0.5*(tilt_A-tilt_B); 
A_spool_t = pi*5e-3*1e-3*1.15;
storke_tilt = 0.6;
V_tt = storke_tilt*tilt_ave*0.9;
at1 = k_tilt/m_tilt;
at2 = b_tilt/m_tilt;
at3 = tilt_ave/m_tilt;
ht1 = beta*(tilt_A^2+tilt_B^2)/V_tt/tilt_ave;
ht2 = beta*(tilt_A^2+tilt_B^2)/V_tt/tilt_ave*C_d*A_spool_t*sqrt(2/rho/(tilt_A^3+tilt_B^3));

u_FFdt = (1/at3*(at1*ydt1+at2*ydt2+ydt3)+ht1*ydt1);
u_FFnt = ht2*sqrt(real(tilt_ave*P_s-sign(u_FFdt).*(tilt_delta*P_s-tilt_ave*(at1*ydt+at2*ydt1+ydt2)/at3)));
u_FFt = u_FFdt./u_FFnt;  % feedforward control input for tilt



%% first case: ILC only 
if option == 1
    k_option_ilc = 1;
    k_option_fb = 0;

        % lift part -----------------------------------------
        xl = zeros(3,tend/dt+1); % state
        xl(:,1) = [0;0;0];  % initial states
        x1l = zeros(3,tend/dt+1); % state
        u_fbl = zeros(1,tend/dt); % feedback control initial values
        ul = zeros(1,tend/dt); % control initial values
        u_ilc_l = zeros(1,tend/dt); % iterative control initial values
        
        k_pl = 10*2; % P gain
        k_il = 0.1*2;  % I gain
        k_ilcl1 = 0.1*0.05/10000;  % iterative control gain 
        k_ilcl2 = 10*4/2*2;  % iterative control gain
        
        err_l = zeros(1,tend/dt+1); % error intitial values
        int_err_l = zeros(1,tend/dt+1); % integral error intitial values
        
        for j = 1:1:n
        
        for i = 1:1:tend/dt
        
        yl(j,i) = xl(1,i);
        err_l(1,i) = ydl(1,i) - yl(j,i);  % tracking error
        
        int_err_l(1,i+1) = err_l(1,i) + int_err_l(1,i);  % integral error
        u_fbl(1,i) = k_pl*err_l(1,i) + k_il*int_err_l(1,i);  % feedback control input
        
        u_ilc_l_j(j,i) =  u_ilc_l(i) + k_ilcl1*(err_l(1,i+1)-err_l(1,i))/dt + k_ilcl2*err_l(1,i); % iterative control input
            
        ul(j,i) = k_option_ilc*u_ilc_l_j(j,i) + k_option_fb*u_fbl(1,i); % ILC + FB
        
        % control input saturation 
        if  ul(j,i) >1 
            ul(j,i) = 1;
        end
        if ul(j,i) <-1 
            ul(j,i) = -1;
        end
        
        xl(:,i+1) = lift_plant(ul(j,i), F_L(i), xl(:,i), dt);
        
        %x1(:,i)= x(:,i); % Store states for next trial
        u_ilc_l(i) = u_ilc_l_j(j,i); % Store input for next trial
        
        end 
        
        x1l(:,i+1)=xl(:,i+1);
        
        end

        
        % tilt part -----------------------------------------
        xt = zeros(3,tend/dt+1); % state
        xt(:,1) = [0;0;0];  % initial states
        x1t = zeros(3,tend/dt+1); % state
        u_fbt = zeros(1,tend/dt); % feedback control initial values
        ut = zeros(1,tend/dt); % control initial values
        u_ilc_t = zeros(1,tend/dt); % iterative control initial values
        
        k_pt = 10*2; % P gain
        k_it = 0.1*2;  % I gain
        k_ilct1 = 0.1*0.05/10000;  % iterative control gain 
        k_ilct2 = 10*4/2*2;  % iterative control gain
        
        err_t = zeros(1,tend/dt+1); % error intitial values
        int_err_t = zeros(1,tend/dt+1); % integral error intitial values
        
        for j = 1:1:n
        
        for i = 1:1:tend/dt
        
        yt(j,i) = xt(1,i);
        err_t(1,i) = ydt(1,i) - yt(j,i);  % tracking error
        
        int_err_t(1,i+1) = err_t(1,i) + int_err_t(1,i);  % integral error
        u_fbt(1,i) = k_pt*err_t(1,i) + k_it*int_err_t(1,i);  % feedback control input
        
        u_ilc_t_j(j,i) =  u_ilc_t(i) + k_ilct1*(err_t(1,i+1)-err_t(1,i))/dt + k_ilct2*err_t(1,i); % iterative control input
            
        ut(j,i) = k_option_ilc*u_ilc_t_j(j,i) + k_option_fb*u_fbt(1,i); % ILC + FB
        
        % control input saturation 
        if  ut(j,i) >1 
            ut(j,i) = 1;
        end
        if ut(j,i) <-1 
            ut(j,i) = -1;
        end
        
        xt(:,i+1) = tilt_plant(ut(j,i), F_T(i), xt(:,i), dt);
        
        %x1(:,i)= xt(:,i); % Store states for next trial
        u_ilc_t(i) = u_ilc_t_j(j,i); % Store input for next trial
        
        end 
        
        x1t(:,i+1)=xt(:,i+1);
        
        end

       
        
end

%% second case: ILC + PI 
if option == 2
    k_option_ilc = 1;
    k_option_fb = 1;

        % lift part -----------------------------------------
        xl = zeros(3,tend/dt+1); % state
        xl(:,1) = [0;0;0];  % initial states
        x1l = zeros(3,tend/dt+1); % state
        u_fbl = zeros(1,tend/dt); % feedback control initial values
        ul = zeros(1,tend/dt); % control initial values
        u_ilc_l = zeros(1,tend/dt); % iterative control initial values
        
        k_pl = 0.01; % P gain
        k_il = 0.1*5*2*2;  % I gain
        k_ilcl1 = 0.1*0.05/10000;  % iterative control gain 
        k_ilcl2 = 10*4/2*2*2;  % iterative control gain
        
        err_l = zeros(1,tend/dt+1); % error intitial values
        int_err_l = zeros(1,tend/dt+1); % integral error intitial values
        
        for j = 1:1:n
        
        for i = 1:1:tend/dt
        
        yl(j,i) = xl(1,i);
        err_l(1,i) = ydl(1,i) - yl(j,i);  % tracking error
        
        int_err_l(1,i+1) = err_l(1,i) + int_err_l(1,i);  % integral error
        u_fbl(1,i) = k_pl*err_l(1,i) + k_il*int_err_l(1,i);  % feedback control input
        
        u_ilc_l_j(j,i) = u_ilc_l(i) + k_ilcl1*(err_l(1,i+1)-err_l(1,i))/dt + k_ilcl2*err_l(1,i); % iterative control input
            
        ul(j,i) = k_option_ilc*u_ilc_l_j(j,i) + k_option_fb*u_fbl(1,i); % ILC + FB
        
        % control input saturation 
        if  ul(j,i) >1 
            ul(j,i) = 1;
        end
        if ul(j,i) <-1 
            ul(j,i) = -1;
        end
        
        xl(:,i+1) = lift_plant(ul(j,i), F_L(i), xl(:,i), dt);
        
        %x1(:,i)= x(:,i); % Store states for next trial
        u_ilc_l(i) = u_ilc_l_j(j,i); % Store input for next trial
        
        end 
        
        x1l(:,i+1)=xl(:,i+1);
        
        end

        
        % tilt part -----------------------------------------
        xt = zeros(3,tend/dt+1); % state
        xt(:,1) = [0;0;0];  % initial states
        x1t = zeros(3,tend/dt+1); % state
        u_fbt = zeros(1,tend/dt); % feedback control initial values
        ut = zeros(1,tend/dt); % control initial values
        u_ilc_t = zeros(1,tend/dt); % iterative control initial values
        
        k_pt = 0.01; % P gain
        k_it = 0.1*5*2*2;  % I gain
        k_ilct1 = 0.1*0.05/10000;  % iterative control gain 
        k_ilct2 = 10*4/2*2*2;  % iterative control gain
        
        err_t = zeros(1,tend/dt+1); % error intitial values
        int_err_t = zeros(1,tend/dt+1); % integral error intitial values
        
        for j = 1:1:n
        
        for i = 1:1:tend/dt
        
        yt(j,i) = xt(1,i);
        err_t(1,i) = ydt(1,i) - yt(j,i);  % tracking error
        
        int_err_t(1,i+1) = err_t(1,i) + int_err_t(1,i);  % integral error
        u_fbt(1,i) = k_pt*err_t(1,i) + k_it*int_err_t(1,i);  % feedback control input
        
        u_ilc_t_j(j,i) =  u_ilc_t(i) + k_ilct1*(err_t(1,i+1)-err_t(1,i))/dt + k_ilct2*err_t(1,i); % iterative control input
            
        ut(j,i) = k_option_ilc*u_ilc_t_j(j,i) + k_option_fb*u_fbt(1,i); % ILC + FB
        
        % control input saturation 
        if  ut(j,i) >1 
            ut(j,i) = 1;
        end
        if ut(j,i) <-1 
            ut(j,i) = -1;
        end
        
        xt(:,i+1) = tilt_plant(ut(j,i), F_T(i), xt(:,i), dt);
        
        %x1(:,i)= xt(:,i); % Store states for next trial
        u_ilc_t(i) = u_ilc_t_j(j,i); % Store input for next trial
        
        end 
        
        x1t(:,i+1)=xt(:,i+1);
        
        end
       
        
end

%% third case: ILC + PI + FB 
if option == 3
    k_option_ilc = 1;
    k_option_fb = 1;
    k_option_ff = 1;

        % lift part -----------------------------------------
        xl = zeros(3,tend/dt+1); % state
        xl(:,1) = [0;0;0];  % initial states
        x1l = zeros(3,tend/dt+1); % state
        u_fbl = zeros(1,tend/dt); % feedback control initial values
        ul = zeros(1,tend/dt); % control initial values
        u_ilc_l = zeros(1,tend/dt); % iterative control initial values
        
        k_pl = 10*2; % P gain
        k_il = 0.1*2;  % I gain
        k_pl = 0.1; % P gain
        k_il = 0.1*5*2*2;  % I gain
        k_ilcl1 = 0.1*0.05/10000;  % iterative control gain 
        k_ilcl2 = 10*4/2*2*2;  % iterative control gain
        
        err_l = zeros(1,tend/dt+1); % error intitial values
        int_err_l = zeros(1,tend/dt+1); % integral error intitial values
        
        for j = 1:1:n
        
        for i = 1:1:tend/dt
        
        yl(j,i) = xl(1,i);
        err_l(1,i) = ydl(1,i) - yl(j,i);  % tracking error
        
        int_err_l(1,i+1) = err_l(1,i) + int_err_l(1,i);  % integral error
        u_fbl(1,i) = k_pl*err_l(1,i) + k_il*int_err_l(1,i);  % feedback control input
        
        u_ilc_l_j(j,i) = u_ilc_l(i) + k_ilcl1*(err_l(1,i+1)-err_l(1,i))/dt + k_ilcl2*err_l(1,i); % iterative control input
            
        ul(j,i) = k_option_ff*u_FFl(1,i) + k_option_ilc*u_ilc_l_j(j,i) + k_option_fb*u_fbl(1,i); % FF + ILC + FB
        
        % control input saturation 
        if  ul(j,i) >1 
            ul(j,i) = 1;
        end
        if ul(j,i) <-1 
            ul(j,i) = -1;
        end
        
        xl(:,i+1) = lift_plant(ul(j,i), F_L(i), xl(:,i), dt);
        
        %x1(:,i)= x(:,i); % Store states for next trial
        u_ilc_l(i) = u_ilc_l_j(j,i); % Store input for next trial
        
        end 
        
        x1l(:,i+1)=xl(:,i+1);
        
        end

        
        % tilt part -----------------------------------------
        xt = zeros(3,tend/dt+1); % state
        xt(:,1) = [0;0;0];  % initial states
        x1t = zeros(3,tend/dt+1); % state
        u_fbt = zeros(1,tend/dt); % feedback control initial values
        ut = zeros(1,tend/dt); % control initial values
        u_ilc_t = zeros(1,tend/dt); % iterative control initial values
        
        k_pt = 0.1; % P gain
        k_it = 0.1*5*2*2;  % I gain
        k_ilct1 = 0.1*0.05/10000;  % iterative control gain 
        k_ilct2 = 10*4/2*2*2;  % iterative control gain
        
        err_t = zeros(1,tend/dt+1); % error intitial values
        int_err_t = zeros(1,tend/dt+1); % integral error intitial values
        
        for j = 1:1:n
        
        for i = 1:1:tend/dt
        
        yt(j,i) = xt(1,i);
        err_t(1,i) = ydt(1,i) - yt(j,i);  % tracking error
        
        int_err_t(1,i+1) = err_t(1,i) + int_err_t(1,i);  % integral error
        u_fbt(1,i) = k_pt*err_t(1,i) + k_it*int_err_t(1,i);  % feedback control input
        
        u_ilc_t_j(j,i) =  u_ilc_t(i) + k_ilct1*(err_t(1,i+1)-err_t(1,i))/dt + k_ilct2*err_t(1,i); % iterative control input
            
        ut(j,i) = k_option_ff*u_FFt(1,i) + k_option_ilc*u_ilc_t_j(j,i) + k_option_fb*u_fbt(1,i); % FF+ ILC + FB
        
        % control input saturation 
        if  ut(j,i) >1 
            ut(j,i) = 1;
        end
        if ut(j,i) <-1 
            ut(j,i) = -1;
        end
        
        xt(:,i+1) = tilt_plant(ut(j,i), F_T(i), xt(:,i), dt);
        
        %x1(:,i)= xt(:,i); % Store states for next trial
        u_ilc_t(i) = u_ilc_t_j(j,i); % Store input for next trial
        
        end 
        
        x1t(:,i+1)=xt(:,i+1);
        
        end
       
        
end


%% plot results
if option == 1
    
    figure()
    subplot(2,1,1)% ---- output
    plot(t(1:end-1), yl(2,:), t(1:end-1), yl(3,:), t(1:end-1), yl(4,:),...
    t(1:end-1), yl(6,:), t(1:end-1), yl(8,:), t(1:end-1), yl(9,:));
    title('ILC lift output for different iteration from 2-9 (mm)');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    subplot(2,1,2) % ---- input
    plot(t(1:end-1), ul(2,:), t(1:end-1), ul(3,:), t(1:end-1), ul(4,:),...
    t(1:end-1), ul(6,:), t(1:end-1), ul(8,:), t(1:end-1), ul(9,:));
    title('ILC lift input for different iterations from 2-9');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    xlabel('time (s)')

    figure()
    subplot(2,1,1) % ---- output iteration 10 vs reference
    plot(t(1:end-1), yl(10,:)*1000, t(1:end-1), ydl(1:end-1)*1000);
    title('Lift refernce v/s output iteration 10 (mm)');
    legend('ILC lift output iteration 10','reference','location','best')
    subplot(2,1,2) % ---- error
    plot(t, err_l*1000);title('ILC lift error [mm]');grid on
    xlabel('time (s)')
        
    figure()
    subplot(2,1,1)% ---- output
    plot(t(1:end-1), yt(2,:), t(1:end-1), yt(3,:), t(1:end-1), yt(4,:),...
    t(1:end-1), yt(6,:), t(1:end-1), yt(8,:), t(1:end-1), yt(9,:));
    title('ILC tilt output for different iteration from 2-9 (mm)');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    subplot(2,1,2) % ---- input
    plot(t(1:end-1), ut(2,:), t(1:end-1), ut(3,:), t(1:end-1), ut(4,:),...
    t(1:end-1), ut(6,:), t(1:end-1), ut(8,:), t(1:end-1), ut(9,:));
    title('ILC tilt input for different iterations from 2-9');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    xlabel('time (s)')

    figure()
    subplot(2,1,1) % ---- output iteration 10 vs reference
    plot(t(1:end-1), yt(10,:)*1000, t(1:end-1), ydt(1:end-1)*1000);
    title('Tilt refernce v/s output iteration 10 (mm)');
    legend('ILC tilt output iteration 10','reference','location','best')
    subplot(2,1,2) % ---- error
    plot(t, err_t*1000);title('ILC tilt error [mm]');grid on
    xlabel('time (s)')
    
    figure()
    subplot(2,1,1)
    plot(t(1:end-1), ul(10,:))
    title('lift ILC input');
    xlabel('time (s)')
    subplot(2,1,2)
    plot(t(1:end-1), ut(10,:))
    title('tilt ILC input');
    legend('total');
    xlabel('time (s)')

end

%% plot results
if option == 2
    
    figure()
    subplot(2,1,1)% ---- output
    plot(t(1:end-1), yl(2,:), t(1:end-1), yl(3,:), t(1:end-1), yl(4,:),...
    t(1:end-1), yl(6,:), t(1:end-1), yl(8,:), t(1:end-1), yl(9,:));
    title('ILC + FB lift output for different iteration from 2-9 (mm)');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    subplot(2,1,2) % ---- input
    plot(t(1:end-1), ul(2,:), t(1:end-1), ul(3,:), t(1:end-1), ul(4,:),...
    t(1:end-1), ul(6,:), t(1:end-1), ul(8,:), t(1:end-1), ul(9,:));
    title('ILC + FB lift input for different iterations from 2-9');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    xlabel('time (s)')

    figure()
    subplot(2,1,1) % ---- output iteration 10 vs reference
    plot(t(1:end-1), yl(10,:)*1000, t(1:end-1), ydl(1:end-1)*1000);
    title('Lift refernce v/s output iteration 10 (mm)');
    legend('ILC + FB lift output iteration 10','reference','location','best')
    subplot(2,1,2) % ---- error
    plot(t, err_l*1000);title('ILC + FB lift error [mm]');grid on
    xlabel('time (s)')
        
    figure()
    subplot(2,1,1)% ---- output
    plot(t(1:end-1), yt(2,:), t(1:end-1), yt(3,:), t(1:end-1), yt(4,:),...
    t(1:end-1), yt(6,:), t(1:end-1), yt(8,:), t(1:end-1), yt(9,:));
    title('ILC + FB tilt output for different iteration from 2-9 (mm)');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    subplot(2,1,2) % ---- input
    plot(t(1:end-1), ut(2,:), t(1:end-1), ut(3,:), t(1:end-1), ut(4,:),...
    t(1:end-1), ut(6,:), t(1:end-1), ut(8,:), t(1:end-1), ut(9,:));
    title('ILC + FB tilt input for different iterations from 2-9');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    xlabel('time (s)')

    figure()
    subplot(2,1,1) % ---- output iteration 10 vs reference
    plot(t(1:end-1), yt(10,:)*1000, t(1:end-1), ydt(1:end-1)*1000);
    title('Tilt refernce v/s output iteration 10 (mm)');
    legend('ILC + FB tilt output iteration 10','reference','location','best')
    subplot(2,1,2) % ---- error
    plot(t, err_t*1000);title('ILC + FB tilt error [mm]');grid on
    xlabel('time (s)')
    
    figure()
    subplot(2,1,1)
    plot(t(1:end-1), ul(10,:),'r',t(1:end-1), u_fbl(1,:),'k',t(1:end-1), u_ilc_l(1,1:end),'b');
    title('lift inputs');
    legend('total','FB','ILC','location','best');
    xlabel('time (s)')
    subplot(2,1,2)
    plot(t(1:end-1), ut(10,:),'r',t(1:end-1), u_fbt(1,:),'k',t(1:end-1), u_ilc_t(1,1:end),'b');
    title('tilt inputs');
    legend('total','FB','ILC','location','best');
    xlabel('time (s)')

end

%% plot results
if option == 3
    
    figure()
    subplot(2,1,1)% ---- output
    plot(t(1:end-1), yl(2,:), t(1:end-1), yl(3,:), t(1:end-1), yl(4,:),...
    t(1:end-1), yl(6,:), t(1:end-1), yl(8,:), t(1:end-1), yl(9,:));
    title('ILC + FF + FB lift output for different iteration from 2-9 (mm)');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    subplot(2,1,2) % ---- input
    plot(t(1:end-1), ul(2,:), t(1:end-1), ul(3,:), t(1:end-1), ul(4,:),...
    t(1:end-1), ul(6,:), t(1:end-1), ul(8,:), t(1:end-1), ul(9,:));
    title('ILC + FF + FB lift input for different iterations from 2-9');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    xlabel('time (s)')

    figure()
    subplot(2,1,1) % ---- output iteration 10 vs reference
    plot(t(1:end-1), yl(10,:)*1000, t(1:end-1), ydl(1:end-1)*1000);
    title('lift refernce v/s output iteration 10 (mm)');
    legend('ILC + FF + FB output iteration 10','reference','location','best')
    subplot(2,1,2) % ---- error
    plot(t, err_l*1000);title('ILC + FF + FB lift error [mm]');grid on
    xlabel('time (s)')
        
    figure()
    subplot(2,1,1)% ---- output
    plot(t(1:end-1), yt(2,:), t(1:end-1), yt(3,:), t(1:end-1), yt(4,:),...
    t(1:end-1), yt(6,:), t(1:end-1), yt(8,:), t(1:end-1), yt(9,:));
    title('ILC + FF + FB tilt output for different iteration from 2-9 (mm)');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    subplot(2,1,2) % ---- input
    plot(t(1:end-1), ut(2,:), t(1:end-1), ut(3,:), t(1:end-1), ut(4,:),...
    t(1:end-1), ut(6,:), t(1:end-1), ut(8,:), t(1:end-1), ut(9,:));
    title('ILC + FF + FB tilt input for different iterations from 2-9');
    legend('Iteratio 2','Iteratio 3','Iteratio 4','Iteratio 6','Iteratio 8','Iteratio 9','location','best');
    xlabel('time (s)')

    figure()
    subplot(2,1,1) % ---- output iteration 10 vs reference
    plot(t(1:end-1), yt(10,:)*1000, t(1:end-1), ydt(1:end-1)*1000);
    title('tilt refernce v/s output iteration 10 (mm)');
    legend('ILC + FF + FB tilt output iteration 10','reference','location','best')
    subplot(2,1,2) % ---- error
    plot(t, err_t*1000);title('ILC + FF + FB tilt error [mm]');grid on
    xlabel('time (s)')
        
    figure()
    subplot(2,1,1)
    plot(t(1:end-1), ul(10,:),'r',t(1:end-1), u_fbl(1,:),'k',t(1:end-1), u_ilc_l(1,1:end),'b',t(1:end-1), u_FFl(1,1:end-1));
    title('lift inputs');
    legend('total','FB','ILC','FF','location','best');
    xlabel('time (s)')
    subplot(2,1,2)
    plot(t(1:end-1), ut(10,:),'r',t(1:end-1), u_fbt(1,:),'k',t(1:end-1), u_ilc_t(1,1:end),'b',t(1:end-1), u_FFt(1,1:end-1));
    title('tilt inputs');
    legend('total','FB','ILC','FF','location','best');
    xlabel('time (s)')

end


%% getting bucket tip position
if option == 1 || option == 2 || option == 3 

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

    % combine lift and tilt trajectory to buckt tip trajectory
    addpath 'G:\My Drive\[Modeling] full order model\modeling with steering control\simulationfiles'
    
    run('BucketDynamicsParameters.m')
    run('SoilParameters.m')
    load('soil_params.mat')
    
    
    % bucket tip trajectory calculation
    [theta_t, x_t, y_t, dot_theta_t, dot_x_t, dot_y_t, ddot_theta_t, ...
        ddot_x_t, ddot_y_t] = ...
        getGeometrywithAcc(ydl+lift_ini*ones(1,length(tnew)), ...
        ydt+tilt_ini*ones(1,length(tnew)), ydl1, ...
        ydt1, ydl2, ydt2, ...
        veh_disp', y_A_init,...
        veh_vel', 0, ...
        veh_acc', 0);
    
    [theta_t_a, x_t_a, y_t_a, dot_theta_t_a, dot_x_t_a, dot_y_t_a, ddot_theta_t_a, ...
        ddot_x_t_a, ddot_y_t_a] = ...
        getGeometrywithAcc(xl(1,1:end)+lift_ini*ones(1,length(tnew)), ...
        xt(1,1:end)+tilt_ini*ones(1,length(tnew)), xl(2,1:end), ...
        xt(2,1:end), ydl2, ydt2, ...
        veh_disp', y_A_init,...
        veh_vel', 0, ...
        veh_acc', 0);

end


%% plot bucket tip trajectory 

if option == 1
    % bucket tip trajectory 
    figure()
    plot(x_t-x_t(1), y_t-y_t(1),'r', "lineWidth", 1.2)
    hold on;
    plot(x_t_a - x_t_a(1), y_t_a - y_t_a(1),'b--',"lineWidth", 1.2)
    legend('Reference','Actual','Location','northwest')
    title("ILC Bucket Tip Trajectory")
    xlabel({'Horizontal Distance (m)'},'FontWeight','bold')
    ylabel({'Vertical Distance (m)'},'FontWeight','bold')

end

if option == 2
    % bucket tip trajectory 
    figure()
    plot(x_t-x_t(1), y_t-y_t(1),'r', "lineWidth", 1.2)
    hold on;
    plot(x_t_a - x_t_a(1), y_t_a - y_t_a(1),'b--',"lineWidth", 1.2)
    legend('Reference','Actual','Location','northwest')
    title("ILC + FB Bucket Tip Trajectory")
    xlabel({'Horizontal Distance (m)'},'FontWeight','bold')
    ylabel({'Vertical Distance (m)'},'FontWeight','bold')

end

if option == 3
    % bucket tip trajectory 
    figure()
    plot(x_t-x_t(1), y_t-y_t(1),'r', "lineWidth", 1.2)
    hold on;
    plot(x_t_a - x_t_a(1), y_t_a - y_t_a(1),'b--',"lineWidth", 1.2)
    legend('Reference','Actual','Location','northwest')
    title("ILC + FF + FB Bucket Tip Trajectory")
    xlabel({'Horizontal Distance (m)'},'FontWeight','bold')
    ylabel({'Vertical Distance (m)'},'FontWeight','bold')

end


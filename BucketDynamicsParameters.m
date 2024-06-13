% clear all;

addpath CommonFuncs 
addpath C:\Users\zhaog\Desktop\BucketLinkage\BucketLinkage
%% Predetermined Values
global x_A_init y_A_init l_AY theta_AY l_AK alpha_K l_AF ...
    gamma_B l_BK l_KF l_AG theta_AG m_Lift W_Lift J_Lift gamma_Lift_cm ...
    l_Lift_cm ...
    beta_E l_DE l_EF gamma_cm l_Lever_cm m_Lever J_Lever W_Lever ...
    l_CD l_BC lambda_Link m_Link J_Link W_Link ...
    alpha_C_i z_t h_t l_Bucket_cm alpha_Bucket_cm J_Bucket m_Bucket ...
    W_Bucket ...
    angle_FAB angle_EDF l_DF;

l_AK = sqrt(1741.60^2 + 25.93^2) * 1e-3;
l_AF = sqrt(1984.15^2 + 823.75^2) * 1e-3;
l_BK = sqrt((2907.14 - 1741.60)^2 + (0.00 - 25.93)^2) * 1e-3;
l_KF = sqrt((1984.15 - 1741.60)^2 + (823.75 - 25.93)^2) * 1e-3;
l_BF = sqrt((2907.14 - 1984.15)^2 + (0.00 - 823.75)^2) * 1e-3;

alpha_K = acos((l_AK^2 + l_AF^2 - l_KF^2) / (2 * l_AK * l_AF));
gamma_B = acos((l_BK^2 + l_KF^2 - l_BF^2) / (2 * l_BK * l_KF));

angle_FAB = atan2(823.75, 1984.15);

l_EF = 0.71120;
l_DF = sqrt(612.46^2 + 563.62^2) * 1e-3;
l_DE = sqrt((612.46 + 711.20)^2 + 563.62^2) * 1e-3;
beta_E = acos((l_DE^2 + l_EF^2 - l_DF^2) / (2 * l_DE * l_EF));

angle_EDF = acos((l_DF^2 + l_DE^2 - l_EF^2) / (2 * l_DE * l_DF));

l_CD = 0.8001;

l_BC = sqrt(0.04334^2 + 0.4258^2);
z_t = 0.03111;
h_t = 1.43599;

x_A_init = (1033.89 - 1010.44) * 1e-3;
y_A_init = (1298.79 - 616.03) * 1e-3;
% x_A_init = -3.763053328257608;
% y_A_init = 1.773590773460247;
l_AY = sqrt(x_A_init^2 + y_A_init^2);
theta_AY = atan(y_A_init / x_A_init);

x_G_init = (1033.89 - 649.62) * 1e-3;
y_G_init = (1298.79 - 616.03) * 1e-3;
l_AG = sqrt((x_G_init - x_A_init)^2 + (y_G_init - y_A_init)^2);
theta_AG = atan((y_G_init - y_A_init) / (x_G_init - x_A_init));

l_Lift_cm = sqrt(1.4894049^2 + 0.24859151^2);
gamma_Lift_cm = atan(248.49 / 1489.56) - atan(25.93 / 1741.60);
l_Lever_cm = sqrt((711.20 - 28.06)^2 + 153.34^2) * 1e-3;
gamma_cm = atan(153.34 / (711.20 - 28.06));
lambda_Link = 0.5;
l_Bucket_cm = sqrt((0.6675)^2 + (0.10358)^2);
alpha_Bucket_cm = atan(425.8 / 43.34) - atan(103.58 / 667.5);

m_Lift = 1572.6209;
m_Lever = 380.2667;
m_Link = 79.144215;
m_Bucket = 1571.1032;
J_Lift = 787.00121;
J_Lever = 77.40479 - m_Lever * (28.063833^2 + 153.33872^2) * 1e-6;
J_Link = 49.697765 - m_Link * 400.05^2 * 1e-6;
J_Bucket = 1.2845283e3 - m_Bucket * l_Bucket_cm^2;
g = 9.81;

W_Lift = m_Lift * g;
W_Lever = m_Lever * g;
W_Link = m_Link * g;
W_Bucket = m_Bucket * g;

l1_ground = 1.472946;  % lift initial value
l2_ground = 1.668526;  % tilt initial value
[xi_B, xi_C] = getBCposition(l1_ground, l2_ground);
alpha_C_i = asin((xi_C - xi_B) / l_BC);

function [x_B, x_C] = getBCposition(q1, q2)
global x_A_init y_A_init l_AY theta_AY l_AK alpha_K l_AF gamma_B l_BK ...
    l_AG theta_AG beta_E l_DE l_EF l_CD l_BC;
    x_G = x_A_init + l_AG * cos(theta_AG);
    y_G = y_A_init + l_AG * sin(theta_AG);
    alpha_A = acos((l_AY^2 + l_AK^2 - q1^2) / (2 * l_AY * l_AK));
    theta_A = alpha_A + theta_AY - pi;
    x_K = x_A_init + l_AK * cos(theta_A);
    y_K = y_A_init + l_AK * sin(theta_A);
    theta_F = theta_A + alpha_K;
    x_F = x_A_init + l_AF * cos(theta_F);
    y_F = y_A_init + l_AF * sin(theta_F);
    gamma_F = atan2(y_F - y_K, x_F - x_K);
    theta_B = gamma_F - gamma_B;
    x_B = x_K + l_BK * cos(theta_B);
    y_B = y_K + l_BK * sin(theta_B);
    l_GF = sqrt((x_G - x_F)^2 + (y_G - y_F)^2);
    alpha_E = acos((l_GF^2 + q2^2 - l_EF^2) / (2 * l_GF * q2));
    theta_F = atan2(y_F - y_G, x_F - x_G);
    theta_E = theta_F + alpha_E;
    x_E = x_G + q2 * cos(theta_E);
    y_E = y_G + q2 * sin(theta_E);
    beta_F = atan2(y_F - y_E, x_F - x_E);
    theta_D = beta_F - beta_E;
    x_D = x_E + l_DE * cos(theta_D);
    y_D = y_E + l_DE * sin(theta_D);
    l_BD = sqrt((x_B - x_D)^2 + (y_B - y_D)^2);
    alpha_D = acos((l_BD^2 + l_CD^2 - l_BC^2) / (2 * l_CD * l_BD));
    theta_B = atan2(y_B - y_D, x_B - x_D);
    theta_C = theta_B + alpha_D;
    x_C = x_D + l_CD * cos(theta_C);
end
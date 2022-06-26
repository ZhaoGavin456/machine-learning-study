function [theta_t, x_t, y_t, dot_theta_t, dot_x_t, dot_y_t, ...
    ddot_theta_t, ddot_x_t, ddot_y_t] = getGeometrywithAcc(q1, q2, ...
    dot_q1, dot_q2, ddot_q1, ddot_q2, x_A, y_A, dot_x_A, dot_y_A, ...
    ddot_x_A, ddot_y_A)
% Calculate the displacement and velocity of the bucket tip

global l_AY theta_AY l_AK alpha_K l_AF gamma_B l_BK l_AG theta_AG ...
    beta_E l_DE l_EF l_CD l_BC alpha_C_i z_t h_t;

alpha_A = acos((l_AY^2 + l_AK^2 - q1.^2) / (2 * l_AY * l_AK));
sin_alpha_A = sin(alpha_A);
cos_alpha_A = cos(alpha_A);
dalpha_Adq1 = q1 ./ (l_AY * l_AK * abs(sin_alpha_A));
d2alpha_Adq1_2 = 1 ./ (l_AY * l_AK * abs(sin_alpha_A)) - ...
    (q1.^2 .* cos_alpha_A) ./ (l_AY^2 * l_AK^2 * abs(sin_alpha_A).^3);

theta_A = alpha_A + theta_AY - pi;
dot_theta_A = dalpha_Adq1 .* dot_q1;
ddot_theta_A = dalpha_Adq1 .* ddot_q1 + d2alpha_Adq1_2 .* dot_q1.^2;

sin_theta_A = sin(theta_A);
cos_theta_A = cos(theta_A);
x_K = x_A + l_AK * cos_theta_A;
y_K = y_A + l_AK * sin_theta_A;
dot_x_K = dot_x_A - l_AK * sin_theta_A .* dot_theta_A;
dot_y_K = dot_y_A + l_AK * cos_theta_A .* dot_theta_A;
ddot_x_K = ddot_x_A - l_AK * (sin_theta_A .* ddot_theta_A + ...
    cos_theta_A .* dot_theta_A.^2);
ddot_y_K = ddot_y_A + l_AK * (cos_theta_A .* ddot_theta_A - ...
    sin_theta_A .* dot_theta_A.^2);

theta_F = theta_A + alpha_K;
dot_theta_F = dot_theta_A;
ddot_theta_F = ddot_theta_A;
sin_theta_F = sin(theta_F);
cos_theta_F = cos(theta_F);
x_F = x_A + l_AF * cos_theta_F;
y_F = y_A + l_AF * sin_theta_F;
dot_x_F = dot_x_A - l_AF * sin_theta_F .* dot_theta_F;
dot_y_F = dot_y_A + l_AF * cos_theta_F .* dot_theta_F;
ddot_x_F = ddot_x_A - l_AF * (sin_theta_F .* ddot_theta_F + ...
    cos_theta_F .* dot_theta_F.^2);
ddot_y_F = ddot_y_A + l_AF * (cos_theta_F .* ddot_theta_F - ...
    sin_theta_F .* dot_theta_F.^2);

gamma_F = atan2(y_F - y_K, x_F - x_K);
theta_B = gamma_F - gamma_B;
l_KF = sqrt((x_K - x_F).^2 + (y_F - y_K).^2);
ptheta_Bpx_F = - (y_F - y_K) ./ l_KF.^2;
ptheta_Bpy_F = (x_F - x_K) ./ l_KF.^2;
p2theta_Bpx_F_2 = 2 * (x_F - x_K) .* (y_F - y_K) ./ l_KF.^4;
p2theta_Bpy_F_2 = - p2theta_Bpx_F_2;
p2theta_Bpx_Fpy_F = ((y_F - y_K).^2 - (x_F - x_K).^2) ./ l_KF.^4;
dot_theta_B = ptheta_Bpx_F .* (dot_x_F - dot_x_K) + ptheta_Bpy_F .* ...
    (dot_y_F - dot_y_K);
ddot_theta_B = ptheta_Bpx_F .* (ddot_x_F - ddot_x_K) + ...
    ptheta_Bpy_F .* (ddot_y_F - ddot_y_K) + ...
    p2theta_Bpx_F_2 .* (dot_x_F - dot_x_K).^2 + ...
    2 * p2theta_Bpx_Fpy_F .* (dot_x_F - dot_x_K) .* (dot_y_F - dot_y_K) ...
    + p2theta_Bpy_F_2 .* (dot_y_F - dot_y_K).^2;

sin_theta_B = sin(theta_B);
cos_theta_B = cos(theta_B);
x_B = x_K + l_BK * cos_theta_B;
y_B = y_K + l_BK * sin_theta_B;
dot_x_B = dot_x_K - l_BK * sin_theta_B .* dot_theta_B;
dot_y_B = dot_y_K + l_BK * cos_theta_B .* dot_theta_B;
ddot_x_B = ddot_x_K - l_BK * (sin_theta_B .* ddot_theta_B + ...
    cos_theta_B .* dot_theta_B.^2);
ddot_y_B = ddot_y_K + l_BK * (cos_theta_B .* ddot_theta_B - ...
    sin_theta_B .* dot_theta_B.^2);

x_G = x_A + l_AG * cos(theta_AG);
y_G = y_A + l_AG * sin(theta_AG);
l_GF = sqrt((x_G - x_F).^2 + (y_G - y_F).^2);
pl_GFpx_G = (x_G - x_F) ./ l_GF;
pl_GFpy_G = (y_G - y_F) ./ l_GF;
p2l_GFpx_G_2 = (y_G - y_F).^2 ./ l_GF.^3;
p2l_GFpx_Gpy_G = - ((x_G - x_F) .* (y_G - y_F)) ./ l_GF.^3;
p2l_GFpy_G_2 = (x_G - x_F).^2 ./ l_GF.^3;
dot_l_GF = pl_GFpx_G .* (dot_x_A - dot_x_F) + pl_GFpy_G .* ...
    (dot_y_A - dot_y_F);
ddot_l_GF = pl_GFpx_G .* (ddot_x_A - ddot_x_F) + pl_GFpy_G .* ...
    (ddot_y_A - ddot_y_F) + p2l_GFpx_G_2 .* (dot_x_A - dot_x_F).^2 + ...
    2 * p2l_GFpx_Gpy_G .* (dot_x_A - dot_x_F) .* (dot_y_A - dot_y_F) + ...
    p2l_GFpy_G_2 .* (dot_y_A - dot_y_F).^2;

alpha_E = acos((l_GF.^2 + q2.^2 - l_EF^2) ./ (2 * l_GF .* q2));
sin_alpha_E = sin(alpha_E);
cos_alpha_E = cos(alpha_E);
palpha_Epl_GF = - 1 ./ abs(sin_alpha_E) .* (1 ./ q2 - cos_alpha_E ./ ...
    l_GF);
palpha_Epq2 = - 1 ./ abs(sin_alpha_E) .* (1 ./ l_GF - cos_alpha_E ./ q2);
p2alpha_Epl_GF_2 = 1 ./ (l_GF .* abs(sin_alpha_E)) .* (1 ./ q2 - ...
    2 * cos_alpha_E ./ l_GF) - cos_alpha_E ./ abs(sin_alpha_E).^3 .* ...
    (1 ./ q2 - cos_alpha_E ./ l_GF).^2;
p2alpha_Epq2_2 = 1 ./ (q2 .* abs(sin_alpha_E)) .* (1 ./ l_GF - ...
    2 * cos_alpha_E ./ q2) - cos_alpha_E ./ abs(sin_alpha_E).^3 .* ...
    (1 ./ l_GF - cos_alpha_E ./ q2).^2;
p2alpha_Epq2pl_GF = 1 ./ abs(sin_alpha_E) .* (1 ./ l_GF.^2 + ...
    1 ./ q2.^2 - cos_alpha_E ./ (q2 .* l_GF)) - cos_alpha_E ./ ...
    abs(sin_alpha_E).^3 .* (1 ./ l_GF - cos_alpha_E ./ q2) .* ...
    (1 ./ q2 - cos_alpha_E ./ l_GF);
dot_alpha_E = palpha_Epq2 .* dot_q2 + palpha_Epl_GF .* dot_l_GF;
ddot_alpha_E = palpha_Epq2 .* ddot_q2 + palpha_Epl_GF .* ddot_l_GF + ...
    p2alpha_Epq2_2 .* dot_q2.^2 + 2 * p2alpha_Epq2pl_GF .* ...
    dot_q2 .* dot_l_GF + p2alpha_Epl_GF_2 .* dot_l_GF.^2;

theta_F = atan2(y_F - y_G, x_F - x_G);
ptheta_Fpx_F = - (y_F - y_G) ./ l_GF.^2;
ptheta_Fpy_F = (x_F - x_G) ./ l_GF.^2;
p2theta_Fpx_F_2 = 2 * (x_F - x_G) .* (y_F - y_G) ./ l_GF.^4;
p2theta_Fpy_F_2 = - p2theta_Fpx_F_2;
p2theta_Fpx_Fpy_F = ((y_F - y_G).^2 - (x_F - x_G).^2) ./ l_GF.^4;
dot_theta_F = ptheta_Fpx_F .* (dot_x_F - dot_x_A) + ptheta_Fpy_F .* ...
    (dot_y_F - dot_y_A);
ddot_theta_F = ptheta_Fpx_F .* (ddot_x_F - ddot_x_A) + ...
    ptheta_Fpy_F .* (ddot_y_F - ddot_y_A) + ...
    p2theta_Fpx_F_2 .* (dot_x_F - dot_x_A).^2 + ...
    2 * p2theta_Fpx_Fpy_F .* (dot_x_F - dot_x_A) .* (dot_y_F - dot_y_A) ...
    + p2theta_Fpy_F_2 .* (dot_y_F - dot_y_A).^2;

theta_E = theta_F + alpha_E;
dot_theta_E = dot_theta_F + dot_alpha_E;
ddot_theta_E = ddot_theta_F + ddot_alpha_E;

sin_theta_E = sin(theta_E);
cos_theta_E = cos(theta_E);
x_E = x_G + q2 .* cos_theta_E;
y_E = y_G + q2 .* sin_theta_E;
dot_x_E = dot_x_A + dot_q2 .* cos_theta_E - q2 .* sin_theta_E .* ...
    dot_theta_E;
dot_y_E = dot_y_A + dot_q2 .* sin_theta_E + q2 .* cos_theta_E .* ...
    dot_theta_E;
ddot_x_E = ddot_x_A + ddot_q2 .* cos_theta_E - 2 * dot_q2 .* ...
    sin_theta_E .* dot_theta_E - q2 .* (sin_theta_E .* ddot_theta_E + ...
    cos_theta_E .* dot_theta_E.^2);
ddot_y_E = ddot_y_A + ddot_q2 .* sin_theta_E + 2 * dot_q2 .* ...
    cos_theta_E .* dot_theta_E + q2 .* (cos_theta_E .* ddot_theta_E - ...
    sin_theta_E .* dot_theta_E.^2);

beta_F = atan2(y_F - y_E, x_F - x_E);
theta_D = beta_F - beta_E;
ptheta_Dpx_F = - (y_F - y_E) / l_EF^2;
ptheta_Dpy_F = (x_F - x_E) / l_EF^2;
p2theta_Dpx_F_2 = 2 * (x_F - x_E) .* (y_F - y_E) ./ l_EF.^4;
p2theta_Dpy_F_2 = - p2theta_Dpx_F_2;
p2theta_Dpx_Fpy_F = ((y_F - y_E).^2 - (x_F - x_E).^2) ./ l_EF.^4;
dot_theta_D = ptheta_Dpx_F .* (dot_x_F - dot_x_E) + ptheta_Dpy_F .* ...
    (dot_y_F - dot_y_E);
ddot_theta_D = ptheta_Dpx_F .* (ddot_x_F - ddot_x_E) + ...
    ptheta_Dpy_F .* (ddot_y_F - ddot_y_E) + ...
    p2theta_Dpx_F_2 .* (dot_x_F - dot_x_E).^2 + ...
    2 * p2theta_Dpx_Fpy_F .* (dot_x_F - dot_x_E) .* (dot_y_F - dot_y_E) ...
    + p2theta_Dpy_F_2 .* (dot_y_F - dot_y_E).^2;

sin_theta_D = sin(theta_D);
cos_theta_D = cos(theta_D);
x_D = x_E + l_DE * cos_theta_D;
y_D = y_E + l_DE * sin_theta_D;
dot_x_D = dot_x_E - l_DE * sin_theta_D .* dot_theta_D;
dot_y_D = dot_y_E + l_DE * cos_theta_D .* dot_theta_D;
ddot_x_D = ddot_x_E - l_DE * (sin_theta_D .* ddot_theta_D + ...
    cos_theta_D .* dot_theta_D.^2);
ddot_y_D = ddot_y_E + l_DE * (cos_theta_D .* ddot_theta_D - ...
    sin_theta_D .* dot_theta_D.^2);

l_BD = sqrt((x_B - x_D).^2 + (y_B - y_D).^2);
pl_BDpx_B = (x_B - x_D) ./ l_BD;
pl_BDpy_B = (y_B - y_D) ./ l_BD;
p2l_BDpx_B_2 = (y_B - y_D).^2 ./ l_BD.^3;
p2l_BDpx_Bpy_B = - (x_B - x_D) .* (y_B - y_D) ./ l_BD.^3;
p2l_BDpy_B_2 = (x_B - x_D).^2 ./ l_BD.^3;
dot_l_BD = pl_BDpx_B .* (dot_x_B - dot_x_D) + pl_BDpy_B .* (dot_y_B - ...
    dot_y_D);
ddot_l_BD = pl_BDpx_B .* (ddot_x_B - ddot_x_D) + ...
    pl_BDpy_B .* (ddot_y_B - ddot_y_D) + ...
    p2l_BDpx_B_2 .* (dot_x_B - dot_x_D).^2 + ...
    2 * p2l_BDpx_Bpy_B .* (dot_x_B - dot_x_D) .* (dot_y_B - dot_y_D) + ...
    p2l_BDpy_B_2 .* (dot_y_B - dot_y_D).^2;

alpha_D = acos((l_BD.^2 + l_CD^2 - l_BC^2) ./ (2 * l_CD * l_BD));
sin_alpha_D = sin(alpha_D);
cos_alpha_D = cos(alpha_D);
palpha_Dpl_BD = - 1 ./ abs(sin_alpha_D) .* (1 / l_CD - cos_alpha_D ./ ...
    l_BD);
p2alpha_Dpl_BD_2 = 1 ./ (l_BD .* abs(sin_alpha_D)) .* (1 / l_CD - 2 * ...
    cos_alpha_D ./ l_BD) - cos_alpha_D ./ abs(sin_alpha_D).^3 .* ...
    (1 / l_CD - cos_alpha_D ./ l_BD).^2;
dot_alpha_D = palpha_Dpl_BD .* dot_l_BD;
ddot_alpha_D = palpha_Dpl_BD .* ddot_l_BD + p2alpha_Dpl_BD_2 .* dot_l_BD.^2;

theta_B = atan2(y_B - y_D, x_B - x_D);
ptheta_Bpx_B = - (y_B - y_D) ./ l_BD.^2;
ptheta_Bpy_B = (x_B - x_D) ./ l_BD.^2;
p2theta_Bpx_B_2 = 2 * (x_B - x_D) .* (y_B - y_D) ./ l_BD.^4;
p2theta_Bpy_B_2 = - p2theta_Bpx_B_2;
p2theta_Bpx_Bpy_B = ((y_B - y_D).^2 - (x_B - x_D).^2) ./ l_BD.^4;
dot_theta_B = ptheta_Bpx_B .* (dot_x_B - dot_x_D) + ptheta_Bpy_B .* ...
    (dot_y_B - dot_y_D);
ddot_theta_B = ptheta_Bpx_B .* (ddot_x_B - ddot_x_D) + ...
    ptheta_Bpy_B .* (ddot_y_B - ddot_y_D) + ...
    p2theta_Bpx_B_2 .* (dot_x_B - dot_x_D).^2 + ...
    2 * p2theta_Bpx_Bpy_B .* (dot_x_B - dot_x_D) .* (dot_y_B - dot_y_D) ...
    + p2theta_Bpy_B_2 .* (dot_y_B - dot_y_D).^2;

theta_C = theta_B + alpha_D;
dot_theta_C = dot_theta_B + dot_alpha_D;
ddot_theta_C = ddot_theta_B + ddot_alpha_D;

sin_theta_C = sin(theta_C);
cos_theta_C = cos(theta_C);
x_C = x_D + l_CD * cos_theta_C;
y_C = y_D + l_CD * sin_theta_C;
dot_x_C = dot_x_D - l_CD * sin_theta_C .* dot_theta_C;
dot_y_C = dot_y_D + l_CD * cos_theta_C .* dot_theta_C;
ddot_x_C = ddot_x_D - l_CD * (sin_theta_C .* ddot_theta_C + ...
    cos_theta_C .* dot_theta_C.^2);
ddot_y_C = ddot_y_D + l_CD * (cos_theta_C .* ddot_theta_C - ...
    sin_theta_C .* dot_theta_C.^2);

alpha_C = atan2(x_C - x_B, y_C - y_B);
palpha_Cpx_C = (y_C - y_B) / l_BC^2;
palpha_Cpy_C = - (x_C - x_B) / l_BC^2;
p2alpha_Cpx_C_2 = - 2 * (y_C - y_B) .* (x_C - x_B) / l_BC^4;
p2alpha_Cpy_C_2 = - p2alpha_Cpx_C_2;
p2alpha_Cpx_Cpy_C = ((x_C - x_B).^2 - (y_C - y_B).^2) / l_BC^4;
dot_alpha_C = palpha_Cpx_C .* (dot_x_C - dot_x_B) + palpha_Cpy_C .* ...
    (dot_y_C - dot_y_B);
ddot_alpha_C = palpha_Cpx_C .* (ddot_x_C - ddot_x_B) + ...
    palpha_Cpy_C .* (ddot_y_C - ddot_y_B) + ...
    p2alpha_Cpx_C_2 .* (dot_x_C - dot_x_B).^2 + ...
    2 * p2alpha_Cpx_Cpy_C .* (dot_x_C - dot_x_B) .* (dot_y_C - dot_y_B) ...
    + p2alpha_Cpy_C_2 .* (dot_y_C - dot_y_B).^2;
theta_t = - alpha_C + alpha_C_i;
dot_theta_t = - dot_alpha_C;
ddot_theta_t = - ddot_alpha_C;

sin_theta_t = sin(theta_t);
cos_theta_t = cos(theta_t);
x_t = x_B + z_t * sin_theta_t + h_t * cos_theta_t;
y_t = y_B - z_t * cos_theta_t + h_t * sin_theta_t;
dot_x_t = dot_x_B + (z_t * cos_theta_t - h_t * sin_theta_t) .* ...
    dot_theta_t;
dot_y_t = dot_y_B + (z_t * sin_theta_t + h_t * cos_theta_t) .* ...
    dot_theta_t;
ddot_x_t = ddot_x_B + (z_t * cos_theta_t - h_t * sin_theta_t) .* ...
    ddot_theta_t - (z_t * sin_theta_t + h_t * cos_theta_t) .* ...
    dot_theta_t.^2;
ddot_y_t = ddot_y_B + (z_t * sin_theta_t + h_t * cos_theta_t) .* ...
    ddot_theta_t + (z_t * cos_theta_t - h_t * sin_theta_t) .* ...
    dot_theta_t.^2;
end


function [F_H, F_V, V_s, d, res2output] = SoilInteraction(bucket_x, ...
    bucket_y, bucket_orient, prev_x, prev_y, prev_d, ...
    V_s_prev, geo_params, soil_params)

global w b g V_heap Lt_max;

alpha = geo_params{1};
ground_x = geo_params{2};
ground_y = geo_params{3};

gamma = soil_params{1};
c = soil_params{2};
c_a = soil_params{3};
phi = soil_params{4};
delta = soil_params{5};
k_c = soil_params{6};
k_phi = soil_params{7};
n = soil_params{8};

if V_s_prev < V_heap
    eta = 1.0;
else
    eta = 0.0;
end

rho = alpha - bucket_orient;
N_g = @(b)(cot(b) - tan(alpha)) * ...
        (cos(alpha) + sin(alpha) * cot(b + phi)) / ...
        (2 * cos(rho + delta) + sin(rho + delta) * cot(b + phi));
beta = fminbnd(N_g, 0, pi/2 - alpha);

dist = sqrt((bucket_x - ground_x)^2 + (bucket_y - ground_y)^2);
theta = atan2((bucket_y - ground_y), (bucket_x - ground_x));
d = dist * sin(alpha - theta);
if d <= 0
    d = 0;
end

d_angle = atan2(bucket_y - prev_y, bucket_x - prev_x);
dh = sqrt((bucket_x - prev_x)^2 + (bucket_y - prev_y)^2) * ...
    cos(alpha - d_angle);

delta_V = dh * (d + prev_d) / 2;
V_s = V_s_prev + eta * delta_V;

if rho < 1e-4
    Lt = Lt_max;
else
    Lt = min(d / sin(rho), Lt_max);
end

denom = cos(rho + delta) + sin(rho + delta) * cot(beta + phi);
N_gamma = (cot(beta) - tan(alpha)) * ...
    (cos(alpha) + sin(alpha) .* cot(beta + phi)) / (2 * denom);
N_c = (1 + cot(beta) * cot(beta + phi)) / denom;
N_a = (sin(rho) - cos(rho) * cot(beta + phi)) / denom;
N_q = (cos(alpha) + sin(alpha) * cot(beta + phi)) / denom;
F = w * (gamma * g * (N_gamma * d^2 + V_s * N_q) + ...
    d * c * N_c + c_a * Lt * N_a);

F_B = w * b * (k_c / b + k_phi) * d^n;
F_H_Bucket = F_B + F * sin(delta) + c_a * w * Lt;
F_V_Bucket = - F * cos(delta);

adhesion = c_a * w * Lt;
gamma_force = w * gamma * g * N_gamma * d^2;
s_force = w * gamma * g * V_s * N_q;
c_force = w * d * c * N_c;
ca_force = w * c_a * Lt * N_a;
cohesion = w * c * d / sin(beta);
w2 = 0.5 * w * gamma * g * d^2 * (cot(beta) - tan(alpha));
wvs = w * gamma * g * V_s;
r = (F * sin(rho + delta) - wvs * sin(alpha) - cohesion * cos(beta) ...
    + adhesion * cos(rho) - w2 * sin(alpha)) / sin(beta + phi);
res2output = {F_H_Bucket, F_V_Bucket, F_B, adhesion, F, ...
    gamma_force, s_force, c_force, ca_force, cohesion, w2, wvs, r};

F_H = F_H_Bucket * cos(bucket_orient) + F_V_Bucket * sin(bucket_orient);
F_V = F_V_Bucket * cos(bucket_orient) - F_H_Bucket * sin(bucket_orient);
% M = F * cos(delta) * L_t / 3;
end


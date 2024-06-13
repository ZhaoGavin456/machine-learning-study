%clear all

global w b g V_max V_heap Lt_max;

w = 2.929;
b = 0.0622;

Lt_max = 1.43599;

g = 9.81;

capacity_stuck = 2.9;
capacity_heaped = 3.4;

V_max = capacity_stuck / w;
V_heap = capacity_heaped / w;

save_fig = false;
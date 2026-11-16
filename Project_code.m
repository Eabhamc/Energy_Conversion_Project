%test commit 2

%__________________________________
%Design Parameters
tube_di = 0.025;
tube_length = 153;

%__________________________________
%Geometry
%__________________________________
%Cavity
c_thickness = 0.1;
c_width = 1.58;
c_length = 0.808;
%__________________________________
% Glass
g_thickness = 0.002;
%__________________________________
%Tube
tube_area = pi * tube_di * tube_length;
%__________________________________
%Plate
num_panels = 12;
Total_p_area = c_width * c_length * num_panels;
inclination_ang = 40/(180*pi);
inclination_ang = 62.43/(180*pi);

%__________________________________
%Water Properties @30 degrees
w_k = 0.615;
w_dyn_vis = 0.798 * 10^-3;
w_Pr = 5.42;
w_cp = 4178;

%__________________________________
%Air Properties @40 degrees
R = 287;
a_density = 1.127;
a_cp = 1007;
a_k = 0.02662;
a_dyn_vis = 1.918*10^-5;
Pr = 0.7255;
atm_p = 10^5;

%__________________________________
%Constants
sigma = 5.67*10^-8;
epsilon_p = 0.9;
alpha_p = 0.9;
accel_grav = 9.81;

%__________________________________



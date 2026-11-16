
%Design Parameters
tube_di = 0.025;
tube_length = 153;
mass_flow_rate = 0.1;
panel_eff = 0.157;

%__________________________________
%Constants
sigma = 5.67*10^-8;
epsilon_p = 0.9;
alpha_p = 0.9;
accel_grav = 9.81;

%__________________________________
%Weather
G = 557.15;
wind_speed = 4.3119;
T_a = 289.672807;
T_sky = 0.0552 * T_a^1.5;
T_w = 283.15;
qs = (alpha_p - panel_eff)* G;
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
total_p_area = c_width * c_length * num_panels;
inclination_ang = 40*(pi/180);
critical_ang = 62.43*(pi/180);

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
a_Pr = 0.7255;
atm_p = 10^5;

%Resistance: Rw, Ra, Rg
%Plate/Tube to water convection
Re_D = (4*mass_flow_rate)/(pi*w_dyn_vis*tube_di);
friction = (0.79*log(Re_D)-1.64)^-2;
if Re_D < 2300
  Nu_D = 3.66;
elseif Re_D >= 3000
  Nu_D = ((friction/8)*(Re_D-1000)*w_Pr)/(1+12.7*(sqrt(friction/8))*(w_Pr^(2/3)-1));
 else
 error('Transition Region');
end
h_w = (Nu_D*w_k)/tube_di;
R_w = 1/(tube_area* h_w);
%Convection from outer glass to outside
h_a =5.7 + 3.8*wind_speed;
R_a = 1/(h_a*total_p_area);
%Glass Conduction
g_k = 0.78;
R_g = g_thickness/(g_k*total_p_area);



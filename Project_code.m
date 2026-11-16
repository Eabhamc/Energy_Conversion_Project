
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

old_Tp = 350;
old_Tg = 300;

for iteration = 1:100
 avg_c_temp = (old_Tp+old_Tg)/2;
 beta = 1/avg_c_temp;
 Ra_c = accel_grav*beta*abs(old_Tp-old_Tg)*c_thickness^3/(a_dyn_vis/(atm_p/(R*avg_c_temp)))^2*a_Pr;


 if Ra_c >= 10^8
   error('outside range');

 elseif Ra_c <= 1708
   Nu_horizontal = 1;

 else
   Nu_horizontal = 1 + 1.44*max(0,(1-1708/Ra_c))+ max(0,((Ra_c^(1/3)/18)-1));

 end

  Nu_vertical = 0.22*(a_Pr*Ra_c/(0.2+a_Pr))^0.28*(c_length/c_thickness)^(-1/4);

  if old_Tp > old_Tg
    Nu_c = Nu_horizontal*(Nu_vertical/Nu_horizontal)^(inclination_ang/critical_ang)*sin(critical_ang)^(inclination_ang/(4*critical_ang));

  elseif old_Tp < old_Tg
    Nu_c = 1+(Nu_vertical-1)*sin(inclination_ang);

  else
    Nu_c = 1;

end


h_c = (Nu_c*a_k)/c_thickness;
R_c = 1/(h_c*total_p_area);
R_top = R_c + R_g + R_a;
R_sky = 1/(epsilon_p*sigma*total_p_area*(old_Tp+T_sky)*(old_Tp^2+T_sky^2));
w_conductance = mass_flow_rate*w_cp*(1-exp(-h_w*tube_area/(mass_flow_rate*w_cp)));

new_Tp = (T_w*w_conductance + T_a/R_top +T_sky/R_sky + qs*total_p_area)/(w_conductance + 1/R_top + 1/R_sky);
new_Tg = new_Tp - (R_c/R_top)*(new_Tp - T_a);

plate_change = abs(old_Tp-new_Tp);
glass_change = abs(old_Tg-new_Tg);

old_Tg = new_Tg;
old_Tp = new_Tp;

if plate_change < 0.001 && glass_change < 0.001
break;
end
end



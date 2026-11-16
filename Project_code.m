
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
weather = [
 0     0          3.244579  285.483667
 1     0          3.190982  285.039667;
 2     0          3.158175  284.702649;
 3     0          3.136035  284.420456;
 4     0          3.101807  284.235246;
 5     0          3.082982  284.084175;
 6    22.804737   3.124526  284.196263;
 7    61.219965   3.233263  284.886597;
 8   162.710526   3.484596  285.875228;
 9   274.003772   3.752456  286.848930;
10   386.304579   3.916684  287.769912;
11   445.779544   4.022912  288.565790;
12   530.500140   4.173035  289.210492;
13   557.153842   4.311947  289.672807;
14   541.482281   4.400368  290.023754;
15   524.346316   4.451368  290.259491;
16   459.854632   4.432719  290.294333;
17   374.063667   4.368035  290.140526;
18   252.689684   4.283105  289.837790;
19   128.558333   4.126860  289.389544;
20    46.958737   3.855632  288.747439;
21    13.445070   3.534579  287.879140;
22     0          3.338491  286.803088;
23     0          3.311509  285.986316;

];


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


%Glass Conduction
g_k = 0.78;
R_g = g_thickness/(g_k*total_p_area);


%Tank
tank_mass = 300;
T_w = 283.15;
end_temp = 45+273.15;
dt = 60;
T_0 = 0;
Max_s = 24*3600;
N_inc = (Max_s/dt);
time = linspace(T_0,Max_s,N_inc+1);
total_seconds = 0;
tank_temp = T_w - 273.15;

old_Tp = 350;
old_Tg = 300;

for jj=1:N_inc
  current_time = time(jj);
  hour = floor(current_time/3600);
  row = hour + 1;

  G = weather(row,2);
  wind_speed = weather(row,3);
  T_a = weather(row, 4);
  T_sky = 0.0552*T_a^1.5;
  qs = (alpha_p - panel_eff)* G;
  %Convection from outer glass to outside
  h_a =5.7 + 3.8*wind_speed;
  R_a = 1/(h_a*total_p_area);
  P_PV = panel_eff*G*total_p_area;

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

T_out = old_Tp + (T_w - old_Tp)*exp(-h_w*tube_area/(mass_flow_rate*w_cp));
Q_useful = mass_flow_rate*w_cp *(T_out-T_w);

if Q_useful < 0
  Q_useful = 0;
end

Q_tank = Q_useful + P_PV;

T_w = T_w+Q_tank*dt/(tank_mass*w_cp);
total_seconds(jj+1) = time(jj+1);
tank_temp(jj+1) = T_w-273.15;

if T_w > end_temp
  break;
end
end

figure (1);
plot(total_seconds/3600, tank_temp);
xlabel('Time of Day [h]');
ylabel('Temperature of the Tank [C]');
grid on;


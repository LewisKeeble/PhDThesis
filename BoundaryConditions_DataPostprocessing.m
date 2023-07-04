%% Chip design
clear

%Design parameters
%Line electrodes have width in x and length in y. Length is equal to the
%sensor_ysize. For no electrodes, best to put electrode_width as zero and halve the
%electrodeisfet_separation
chip.N_x = 64;
chip.N_y = 64; %N_x x N_y array
chip.electrode_width = 0E-6;
chip.electrode_thickness = 0E-6;
chip.N_electrodes = 0;%chip.N_x+1; %Number of electrodes on chip

chip.electrodeisfet_separation = 0E-6; %Separation between the outermost point of electrodes and the edges of ISFETs along x axis
chip.isfetisfet_separation = 2E-6; %Separation between the edges of neighbouring ISFETs in the y axis
chip.isfet_width = 9.5E-6;
chip.isfet_length = 10.5E-6;

%Distance between first and last ISFETs and edge of chip in y-axis
chip.sensor_startSeparation = 0E-6;
chip.sensor_endSeparation = 0E-6;

%Width of individual castellations in the y-axis. Extension in x-axis
chip.castellation_width = 0;
chip.castellation_extension = 0;

%If only 1 castellation per pixel, castellation separation will be set to 0
chip.ncastellation_perpixel = 0;
chip.castellation_separation = 0;

%Separation between the edge of the chip and the wall of the reaction
%chamber
chip.wall_separation_xpos = (1.75e-3 + 2*0.262e-3 - (chip.N_x*(2*chip.electrodeisfet_separation+chip.isfet_width) + chip.N_electrodes*chip.electrode_width))/2;
chip.wall_separation_xneg = chip.wall_separation_xpos;
chip.wall_separation_ypos = 0;
chip.wall_separation_yneg = 1.5e-3 - (chip.sensor_startSeparation + chip.N_y*(chip.isfet_length) + (chip.N_y-1)*(chip.isfetisfet_separation) + chip.sensor_endSeparation);

%Calculate the total x and y size of the sensor array
sensor_xsize = chip.N_x*(2*chip.electrodeisfet_separation+chip.isfet_width) + chip.N_electrodes*chip.electrode_width;
sensor_ysize = chip.sensor_startSeparation + chip.N_y*(chip.isfet_length) + (chip.N_y-1)*(chip.isfetisfet_separation) + chip.sensor_endSeparation;

reaction_xsize = sensor_xsize + chip.wall_separation_xpos + chip.wall_separation_xneg;
reaction_ysize = sensor_ysize + chip.wall_separation_ypos + chip.wall_separation_yneg;

solution_vol = 12*10^(-9); %12uL
solution_height = 2.5e-3;
sensor_vol = sensor_xsize*sensor_ysize*solution_height;

%% Load Data
sensor_nH_abs = load('SensornH_Kenny1e3_BurstFactor1e2_Absorbing_2.txt');
sensor_nH_ref = load('SensornH_Kenny1e3_BurstFactor1e2_Reflecting_2.txt');
protons_released = load('ProtonsReleased_Kenny1e3_BurstFactor1e2_2.txt');

%% Postprocess Data

% Constants
threshold.temperature = 336.15; %In Kelvin. 336.15K for 63C for LAMP reaction
threshold.k_B = 1.381E-23; %Boltzmann's constant
threshold.e = 1.602E-19; %Electron charge
threshold.permittivity_water = 65; %Roughly 65 at 63C
threshold.permittivity_vacuum = 8.854E-12;
threshold.d_OHP = 1E-6; %Rough estimate
threshold.n_conc = 2.1; %In mol/L (M). Concentration of charged species
C_stern = 20e-6/1e-2^2;

% From 1D COMSOL simulation
c_SiO_i = 1.6562e-7;
c_SiO_f = 1.4741e-8;

attenuation_factor = 0.48;

% pH change
pH_0 = 8.88;
pH_1 = 7.73;

% Target change in # protons
% Includes surface bound protons
BurstFactor = 0.01;
n_H_target = (10^(-pH_1)-10^(-pH_0))*6.022e23*1000*reaction_xsize*reaction_ysize*solution_height + abs(c_SiO_f-c_SiO_i)*reaction_xsize*reaction_ysize*6.022e23;
nH_actual = max(protons_released/BurstFactor,[],2);
output_scale = n_H_target./nH_actual;
output_scale = repmat(output_scale,1,chip.N_x*chip.N_y,2401);
output_scale = permute(output_scale,[2 3 1]);

%{
sigma = ones(chip.N_x*chip.N_y,length(protons_released),size(sensor_nH,3))*-1.602e-19*c_SiO_i*6.022e23;
sigma = sigma - cumsum(sensor_nH,2)*-1.602e-19/(chip.isfet_width*chip.isfet_length*BurstFactor);
V_T = 2*threshold.k_B*threshold.temperature*asinh(-sigma/ ...
(sqrt(8*threshold.k_B*threshold.temperature*threshold.permittivity_water*threshold.permittivity_vacuum*threshold.n_conc))) ...
+ sigma/C_stern;

figure(1)
hold on
for i = 1:size(sensor_nH,3)
plot(mean(V_T(:,:,i)))
end
hold off
shg
%}

% Absorbing
sigma_scaled_abs = ones(chip.N_x*chip.N_y,length(protons_released))*-1.602e-19*c_SiO_i*6.022e23;
sigma_scaled_abs = sigma_scaled_abs - cumsum(sensor_nH_abs.*output_scale/BurstFactor,2)*-1.602e-19/(chip.isfet_width*chip.isfet_length);
V_T_scaled_abs = (2*threshold.k_B*threshold.temperature*asinh(-sigma_scaled_abs/ ...
(sqrt(8*threshold.k_B*threshold.temperature*threshold.permittivity_water*threshold.permittivity_vacuum*threshold.n_conc))) ...
+ sigma_scaled_abs/C_stern)*attenuation_factor;

% Reflecting
sigma_scaled_ref = ones(chip.N_x*chip.N_y,length(protons_released))*-1.602e-19*c_SiO_i*6.022e23;
sigma_scaled_ref = sigma_scaled_ref - cumsum(sensor_nH_ref.*output_scale/BurstFactor,2)*-1.602e-19/(chip.isfet_width*chip.isfet_length);
V_T_scaled_ref = (2*threshold.k_B*threshold.temperature*asinh(-sigma_scaled_ref/ ...
(sqrt(8*threshold.k_B*threshold.temperature*threshold.permittivity_water*threshold.permittivity_vacuum*threshold.n_conc))) ...
+ sigma_scaled_ref/C_stern)*attenuation_factor;

% Experimental
t_Kenny = load('t_Kenny');
V_Kenny = load('V_Kenny');

figure(3)
hold on
plot(t_Kenny.t_Kenny,V_Kenny.V_Kenny*1e-3)
plot(1:2401,mean(V_T_scaled_abs-V_T_scaled_abs(1)))
plot(1:2401,mean(V_T_scaled_ref-V_T_scaled_ref(1)))
hold off
shg
xlabel('Time (s)')
ylabel('Threshold Voltage Change (V)')
legend('Exp','Abs','Ref')
savefig('BoundaryConditions_ExpVsSim')

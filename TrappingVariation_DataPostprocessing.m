
%clear

%% Load Data
path = ["_Titan_Trapping_1e3Copy_BurstFactor1e1","_Titan_Trapping_1e3Copy_BurstFactor1e1","_Titan_Trapping_1e3Copy_BurstFactor1e2","_Titan_Trapping_1e3Copy_BurstFactor1e2","_Titan_Trapping_1e3Copy_BurstFactor1e2"];
dir = ["0","25","50","75","100"];
%{
for j = 5
    disp('protons_released')
    protons_released = [];
    sensor_nH = [];
    run_error = [];
    for i = 1:50
        disp(i)
        try
            input = load(sprintf('%s\\ProtonsReleased%s_Run%d',dir(j),path(j),i));
        catch
            run_error(end+1) = i;
            disp('run error\n')
            continue
        end
        input_size = size(input);
        if input_size(1) == 1 && input_size(2) == 2401 && length(input_size) == 2
            protons_released(end+1,:) = input;
        else
            run_error(end+1) = i;
            disp('run error\n')
        end
    end

    disp('sensor_nH')
    run = 1:50;
    run(run_error)=[];
    for i = run
        disp(i)
        try
            input = load(sprintf('%s\\SensornH%s_Run%d.txt',dir(j),path(j),i));
        catch
            run_error(end+1) = i;
            disp('run error\n')
            protons_released(i,:) = [];
            continue
        end
        input_size = size(input);
        if input_size(1) == 10000 && input_size(2) == 2401 && length(input_size) == 2 && isempty(find(run_error==i,1))
            if i == 1
                sensor_nH = input;
            else
                sensor_nH(:,:,end+1) = input;
            end
        else
            run_error(end+1) = i;
            disp('run error\n')
            protons_released(i,:) = [];
        end
    end
    save(sprintf('%s\\SensornH%s',dir(j),path(j)),'sensor_nH','-v7.3')
    save(sprintf('%s\\ProtonsReleased%s',dir(j),path(j)),'protons_released','-v7.3')
    clear('sensor_nH')
    clear('protons_released')
end
%}
%% Postprocess Data
%% Chip design
%Design parameters
%Line electrodes have width in x and length in y. Length is equal to the
%sensor_ysize. For no electrodes, best to put electrode_width as zero and halve the
%electrodeisfet_separation
chip.N_x = 100;
chip.N_y = 100; %N_x x N_y array
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

output_scale_final = 4.814076274732945;

% Target change in # protons
BurstFactor = 0.01;
n_H_target = (10^(-pH_1)-10^(-pH_0))*6.022e23*1000*reaction_xsize*reaction_ysize*solution_height + abs(c_SiO_f-c_SiO_i)*reaction_xsize*reaction_ysize*6.022e23;
%{
for i = 1:5
    %{
    protons_released = load(sprintf('%s\\ProtonsReleased%s.mat',dir(i),path(i)));
    nH_actual = max(protons_released.protons_released/BurstFactor,[],2);
    output_scale = permute(repmat(n_H_target./nH_actual,1,chip.N_x*chip.N_y,size(protons_released.protons_released,2)),[2 3 1]);
    clear('nH_actual');
    clear('protons_released');
    %}
    %{
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
    sensor_nH = load(sprintf('%s\\SensornH%s.mat',dir(i),path(i)));
    if sum(sum(sensor_nH.sensor_nH(:,:,1)))==0
        sensor_nH(:,:,1)=[];
    end
    sigma_scaled = ones(chip.N_x*chip.N_y,size(sensor_nH.sensor_nH,2),size(sensor_nH.sensor_nH,3))*-1.602e-19*c_SiO_i*6.022e23;
    sigma_scaled = sigma_scaled - cumsum(sensor_nH.sensor_nH/BurstFactor,2)*-1.602e-19/(chip.isfet_width*chip.isfet_length);
    clear('output_scale');
    clear('sensor_nH');
    V_T_scaled = (2*threshold.k_B*threshold.temperature*asinh(-sigma_scaled/ ...
    (sqrt(8*threshold.k_B*threshold.temperature*threshold.permittivity_water*threshold.permittivity_vacuum*threshold.n_conc))) ...
    + sigma_scaled/C_stern)*attenuation_factor;
    save(sprintf('%s\\V_T_unscaled%s',dir(i),path(i)),'V_T_scaled','-v7.3')
%}
end
%}

V_T = cell(1,5);
path = ["_Titan_Trapping_1e3Copy_BurstFactor1e1","_Titan_Trapping_1e3Copy_BurstFactor1e1","_Titan_Trapping_1e3Copy_BurstFactor1e2","_Titan_Trapping_1e3Copy_BurstFactor1e2","_Titan_Trapping_1e3Copy_BurstFactor1e2"];
dir = ["0","25","50","75","100"];

for i = 1:5
input = load(sprintf('%s\\V_T_unscaled%s.mat',dir(i),path(i)));
V_T_temp = input.V_T_scaled;
onpixels = find(V_T_temp(:,end)-V_T_temp(:,1)>0);
V_T{i} = squeeze(mean(V_T_temp(onpixels,:,:)));
end
clear('input')
%}
plot_colour = ['-g';'-b';'-r';'-m';'-c'];
figure(1)
hold on
for i=[1,5]
   shadedErrorBar(0:length(V_T{i})-1,(V_T{i}-V_T{i}(1,:))',{@mean,@std},'lineprops',plot_colour(i)) 
end
legend('0%','25%','50%','75%','100%')
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Threshold Voltage (V)')
shg

axes('Position',[.6 .2 .3 .3])
box on
for i=[1,5]
   shadedErrorBar(0:length(V_T{i})-1,(V_T{i}-V_T{i}(1,:))',{@mean,@std},'lineprops',plot_colour(i)) 
end
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
xlim([300 600])
box off
hold off
shg

%[h_m h_i]=inset(h_wide,h_narrow);

t=0:2400;
ttp = cell(1,5);
for i = 1:5
    ttp{i} = zeros(1,size(V_T{i},2));
    for j = 1:size(V_T{i},2)
        % Plot TTP as the time of the peak derivative
        [~,max_diff_idx] = max(movmean(diff(V_T{i}(:,j)-V_T{i}(1,j)),50));

        % Plot TTP as the time when of the amplification line crossing y=0
        ind_range = 5;
        xData=t(max_diff_idx-ind_range:max_diff_idx+ind_range)/60;
        yData=movmean(V_T{i}(:,j)-V_T{i}(1,j),50);
        yData=yData(max_diff_idx-ind_range:max_diff_idx+ind_range);
        coeffs=coeffvalues(fit(xData',yData,'poly1'));
        amplif_fit_x = 7:22;
        amplif_fit=coeffs(1)*amplif_fit_x+coeffs(2);
        ttp{i}(j) = -coeffs(2)/coeffs(1);
    end
end

ttp_av = zeros(1,5);
ttp_std = zeros(1,5);
for i = 1:5
    ttp_av(i) = mean(ttp{i});
    ttp_std(i) = std(ttp{i});
end

p_fit = polyfit([0,25,50,75,100],ttp_av(1:5),1);
y_fit = polyval(p_fit,[0,25,50,75,100]);
ss_tot = sum((ttp_av(1:5)-mean(ttp_av)).^2);
ss_res = sum((ttp_av(1:5)-y_fit).^2);
r_sq = 1-ss_res/ss_tot;

figure(4)
hold on
plot([0,25,50,75,100],y_fit,'-b')
errorbar([0,25,50,75,100],ttp_av(1:5),ttp_std(1:5),"o","Color","blue",'LineWidth',2)
hold off
legend(sprintf('ISFET -\ny=-%0.2fx+%0.2f\nr^2=%0.4f',p_fit(1),p_fit(2),r_sq))
xlabel('log_1_0(DNA Copy)')
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
ylabel('ttp (mins)')
shg


%{
figure(3)
hold on
for i = 1:5
   plot(0:length(V_T{i})-2,diff(mean(V_T{i}-V_T{i}(1,:),2) ))
end
legend('0','25','50','75','100')
xlabel('Time (s)')
ylabel('Threshold Voltage (V)')
hold off
shg
%}
clear

channel = 2;

start_idx = 1267;
end_idx = 1900;
offset_idx = 32;

Eoff_start_idx = 290;
Eoff_end_idx = 383;

Eon_line_start_idx = 866;
Eon_line_end_idx = 1559;

x = 0:210e-6/(end_idx-start_idx):210e-6;

Eon = importdata('Titan4_400pMLambda_2V_500Hz_EOn_3mins.csv');
Eon = Eon.data(2+start_idx:end_idx,channel+3);
Eon = Eon./Eon(1);
%Eon = Eon-Eon(1);
Eon_av = max(Eon)*ones(1,length(Eon));

Eoff = importdata('Titan4_400pMLambda_2V_500Hz_EOff.csv');
Eoff = Eoff.data(2+start_idx-offset_idx:end_idx-offset_idx,channel+3);
Eoff = Eoff./Eoff(1);
%Eoff = Eoff/max(Eon);
%Eoff = Eoff-Eoff(1);
Eoff_av = mean(Eoff(Eoff_start_idx:Eoff_end_idx))*ones(1,length(Eon));

Eon_line = importdata('Titan3_400pMLambda_2V_500Hz_EOn_3mins.csv');
Eon_line = Eon_line.data(2+Eon_line_start_idx+offset_idx:Eon_line_end_idx+offset_idx,channel+3);
%Eon_line = Eon_line-Eon_line(3);
Eon_line = Eon_line./Eon_line(3);
Eon_line_av = max(Eon_line)*ones(1,length(Eon));

%Eon = Eon/max(Eon);



figure(1)
clf
hold on
plot(x(1:length(Eon)),Eoff,'LineWidth',2)
plot(x(1:length(Eon)),Eon,'LineWidth',2)
plot(x(1:length(Eon)),Eon_line(1:length(Eon)),'LineWidth',2)
plot(x(1:length(Eon)), Eoff_av, '--k')
plot(x(1:length(Eon)), Eon_line_av, '--k')
plot(x(1:length(Eon)), Eon_av, '--k')
hold off
set(gca,'FontSize',10)
set(gca,'LineWidth',1)
xlim([0 2e-4])
ylim([1 2.4])
xlabel('Distance (m)')
ylabel('Normalised Fluorescence')
legend('Control','Sawtooth','Line')
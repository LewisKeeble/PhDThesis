clear
input_chip = importdata('TAHU4_2_Polarisation_2_5uA_50nAScanRate_2.txt');
%input_chip = importdata('TAHU4_2_Polarisation_2_5uA_50nAScanRate_1.txt');
%input = importdata('TAHU4_2_Polarisation_20uA_100nAScanRate.txt');
input_slide = importdata('AgClSlide_100mA_1mAScanRate.txt');
input_wire = importdata('AgClWire_Polarisation_3mA_10uAScanRate.txt');
data_chip = input_chip.data;
data_slide = input_slide.data;
data_wire = input_wire.data;

electrode_area_chip = (90e-6*90e-6+4*90e-6*15e-6)*100; %cm^2
electrode_area_slide = 25e-3*10e-3*100; %cm^2
electrode_area_wire = 25e-3*pi*1e-3*100; %cm^2

figure(1)
clf
plot(abs(data_chip(:,2))/electrode_area_chip,data_chip(:,1),'LineWidth',1)
set(gca, 'XScale', 'log')
set(gca,'FontSize',10)
xlabel('Current (A/cm^2)')
ylabel('Potential (V)')

figure(2)
clf
hold on
plot(abs(data_chip(:,2))/electrode_area_chip,data_chip(:,1),'LineWidth',1)
plot(abs(data_slide(:,2))/electrode_area_slide,data_slide(:,1),'LineWidth',1)
plot(abs(data_wire(:,2))/electrode_area_wire,data_wire(:,1),'LineWidth',1)
set(gca, 'XScale', 'log')
set(gca,'FontSize',10)
xlabel('Current (A/cm^2)')
ylabel('Potential (V)')
legend('Pad','Slide','Wire')
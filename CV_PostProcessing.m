clear

%% Ag/AgCl glass vs plated RE
%{
input_commercial = importdata('TAHU2_1_AuISEWE_PtCE_AgAgClRE_25mVSweep2_1mMK4FeC6N6_100uMK3FeC6N6_100mMKCl.txt');
input_electroplated = importdata('TAHU2_1_AuISEWE_PtCE_ElectroplatedRE_25mVSweep_1mMK4FeC6N6_100uMK3FeC6N6_100mMKCl.txt');

figure(1)
hold on
plot(input_commercial.data(:,1),input_commercial.data(:,2))
plot(input_electroplated.data(:,1),input_electroplated.data(:,2))
hold off
legend('Commercial glass RE','Electroplated RE')
xlabel('Voltage (V)')
ylabel('Current (A)')
set(gca,'FontSize',12)
set(gca,'LineWidth',1)
%}
electrode_area_chip_AgCl = (90e-6*90e-6+4*90e-6*15e-6)*10000; %cm^2
electrode_area_chip_ENIG = (80e-6*80e-6+4*90e-6*5e-6)*10000; %cm^2
electrode_area_slide = 25e-3*10e-3*10000; %cm^2
electrode_area = electrode_area_chip_AgCl;
%% Peak current vs scan rate

scanrate = [10,20,30,40,50,60,70,80,90,100];
input_CV = {};
for i = 1:length(scanrate) 
   input_CV_temp = importdata(sprintf('TAHU4_2_AgCl_ISE_Increasing//%dmVSweep.txt', scanrate(i))); 
   %input_CV_temp = importdata(sprintf('AgClSlide_CV_Increasing//%dmVSweep.txt', scanrate(i)));
   %input_CV_temp = importdata(sprintf('AgClWire_CV_Increasing//%dmVSweep.txt', scanrate(i)));
   %input_CV_temp = importdata(sprintf('TAHU4_2_ENIG_ISE_Increasing//%dmVScan.txt', scanrate(i))); 
   %input_CV_temp = importdata(sprintf('TAHU4_1_WE_CV_Increase//%dmVScan.txt', scanrate(i))); 
   %input_CV_temp = importdata(sprintf('GoldSlide_CV_Increasing//GoldSlide_CV_%dmVSweep.txt', scanrate(i))); 
   input_CV{end+1} = input_CV_temp.data;
end

figure(2)
hold on
for i = 1:length(scanrate)
   %input_CV{i}(:,2)=input_CV{i}(:,2)-input_CV{i}(1,2);
   plot(input_CV{i}(:,1),input_CV{i}(:,2)./electrode_area,'LineWidth',2) 
end
xlabel('Voltage (V)')
ylabel('Current (A)')
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
%legend(string(scanrate))
hold off

figure(5)
hold on
for i = 1:length(scanrate)
   colour = 0.8 - 0.8*i/(length(scanrate));
   plot(input_CV{i}(:,1),movmean(input_CV{i}(:,2),5)./electrode_area,'Linewidth',2,'Color',[1 colour colour]) 
end
xlabel('Voltage (V)')
ylabel('Current Density (A/cm^2)')
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
xlim([0 0.5])
%legend(string(scanrate))
hold off

% Cathodic Peaks 

cathodicpeaks_val = zeros(1,length(scanrate));
cathodicpeaks_idx = zeros(1,length(scanrate));
for i = 1:length(scanrate)
   [~,switchover_idx] = max(input_CV{i}(:,1)); 
   [cathodicpeaks_val(i),temp_idx] = min(input_CV{i}(switchover_idx:end,2)); 
   cathodicpeaks_idx(i) = temp_idx + switchover_idx;
end

cathodicpeaks_fit_1 = polyfit(sqrt(scanrate),cathodicpeaks_val,1);
y_fit_c = polyval(cathodicpeaks_fit_1,sqrt(scanrate));
SStot = sum((cathodicpeaks_val-mean(cathodicpeaks_val)).^2);                    % Total Sum-Of-Squares
SSres_1 = sum((cathodicpeaks_val-y_fit_c).^2);                       % Residual Sum-Of-Squares
Rsq_c = 1-SSres_1/SStot;

cathodicpeaks_fit_2 = polyfit(sqrt(scanrate),cathodicpeaks_val,2);
y_fit_2 = polyval(cathodicpeaks_fit_2,sqrt(scanrate));
SSres_2 = sum((cathodicpeaks_val-y_fit_2).^2);                       % Residual Sum-Of-Squares
Rsq_2 = 1-SSres_2/SStot;

figure(3)
hold on 
plot(sqrt(scanrate),cathodicpeaks_val./electrode_area,'LineWidth',1)
plot(sqrt(scanrate),y_fit_c./electrode_area,'LineWidth',1)
%plot(sqrt(scanrate),y_fit_2)
hold off
xlabel('scanrates^1^/^2 (mV/s^1^/^2)')
ylabel('Cathodic Peak Current (A)')
legend('Data',sprintf('Linear Fit - r^2=%0.4f',Rsq_c),sprintf('Square Fit - r^2=%0.4f',Rsq_2))
set(gca,'FontSize',10)

%Anodic Peaks

scanrate = [10,20,30,40,50,60,70,80,90,100];
anodicpeaks_val = zeros(1,length(scanrate));
anodicpeaks_idx = zeros(1,length(scanrate));
for i = 1:length(scanrate)
   [anodicpeaks_val(i),anodicpeaks_idx(i)] = max(input_CV{i}(:,2)); 
end

anodicpeaks_fit_1 = polyfit(sqrt(scanrate),anodicpeaks_val,1);
y_fit_a = polyval(anodicpeaks_fit_1,sqrt(scanrate));
SStot = sum((anodicpeaks_val-mean(anodicpeaks_val)).^2);                    % Total Sum-Of-Squares
SSres_1 = sum((anodicpeaks_val-y_fit_a).^2);                       % Residual Sum-Of-Squares
Rsq_a = 1-SSres_1/SStot;

anodicpeaks_fit_2 = polyfit(sqrt(scanrate),anodicpeaks_val,2);
y_fit_2 = polyval(anodicpeaks_fit_2,sqrt(scanrate));
SSres_2 = sum((anodicpeaks_val-y_fit_2).^2);                       % Residual Sum-Of-Squares
Rsq_2 = 1-SSres_2/SStot;

figure(4)
hold on 
plot(sqrt(scanrate),anodicpeaks_val./electrode_area,'LineWidth',1)
plot(sqrt(scanrate),y_fit_a./electrode_area,'LineWidth',1)
%plot(sqrt(scanrate),y_fit_2)
hold off
xlabel('scanrates^1^/^2 (mV/s^1^/^2)')
ylabel('Cathodic Peak Current (A)')
legend('Data',sprintf('Linear Fit - r^2=%0.4f',Rsq_a),sprintf('Square Fit - r^2=%d',Rsq_2))
set(gca,'FontSize',10)

figure(6)
hold on 
plot(sqrt(scanrate),anodicpeaks_val./electrode_area,'-b','LineWidth',2)
plot(sqrt(scanrate),abs(cathodicpeaks_val)./electrode_area,'-r','LineWidth',2)
plot(sqrt(scanrate),y_fit_a./electrode_area,'--b','LineWidth',2)
plot(sqrt(scanrate),abs(y_fit_c)./electrode_area,'--r','LineWidth',2)
hold off
xlabel('scanrate^1^/^2 (mV/s^1^/^2)')
ylabel('Peak Current Density (A/cm^2)')
legend(sprintf('Anodic Peaks - r^2=%0.4f',Rsq_a),sprintf('Cathodic Peaks - r^2=%0.4f',Rsq_c))
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
xlim([3 10])

figure(7)
hold on
for i = 1:length(scanrate)
   plot(input_CV{i}(:,1),input_CV{i}(:,2)) 
   plot(input_CV{i}(anodicpeaks_idx(i),1),input_CV{i}(anodicpeaks_idx(i),2),'ok')
   plot(input_CV{i}(cathodicpeaks_idx(i),1),input_CV{i}(cathodicpeaks_idx(i),2),'ok')
end
xlabel('Voltage (V)')
ylabel('Current (A)')
%legend(string(scanrate))
hold off


%Peak-to-peak
peaktopeak = zeros(1,length(scanrate));
for i = 1:length(scanrate)
    peaktopeak(i) = input_CV{i}(anodicpeaks_idx(i),1) - input_CV{i}(cathodicpeaks_idx(i),1);
end
mean(peaktopeak)

%Peak ratio
i_lambda = zeros(1,length(scanrate));
for i = 1:length(scanrate)
    [~,switch_idx] = max(input_CV{i}(:,1));
    i_lambda(i) = input_CV{i}(switch_idx,2);
end
peakratio = abs(cathodicpeaks_val)./anodicpeaks_val + 0.48*i_lambda./anodicpeaks_val + 0.086;
mean(peakratio)

%FWHM
halfmax = anodicpeaks_val/2;
fwhm = zeros(1,length(scanrate));
for i = 1:length(scanrate)
    idx_1 = find(input_CV{i}(:,2)>=halfmax(i),1,'first');
    idx_2 = find(input_CV{i}(:,2)>=halfmax(i),1,'last');
    fwhm(i) = input_CV{i}(idx_2,1)-input_CV{i}(idx_1,1);
end
mean(fwhm)

figure(8)
input_temp = importdata('TAHU2_1_AuISEWE_AuPlatedCE_ElectroplatedRE_25mVSweep_1mMK4FeC6N6_100uMK3FeC6N6_100mMKCl.txt'); 
input_CV_onchip = input_temp.data;
input_temp = importdata('TAHU2_1_AuISEWE_PtCE_AgAgClRE_25mVSweep3_1mMK4FeC6N6_100uMK3FeC6N6_100mMKCl.txt'); 
input_CV_offchip = input_temp.data;
hold on
plot(input_CV_onchip(:,1),movmean(input_CV_onchip(:,2),10),'LineWidth',1) 
plot(input_CV_offchip(:,1),movmean(input_CV_offchip(:,2),10),'LineWidth',1) 
xlabel('Voltage (V)')
ylabel('Current (A)')
legend('on-chip','external')
set(gca,'FontSize',10)
xlim([0 0.3])
hold off

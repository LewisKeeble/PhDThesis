clear
input = importdata('TAHU4_2_OCP.txt');
%input = importdata('TAHU2_1_OCP_100mMKCl.txt');
%input = importdata('TAHU2_3_OCP_500mMKCL_7_5minsFeCl3.txt');
%input = importdata('TAHU_OCP_500mMKCl.txt');
%input = importdata('TAHU_OCP_500mMKCl_10min50mMFeCl3.txt');
%input = importdata('TAHU_OCP_500mMKCl_10min50mMFeCl3.txt');
OCP_raw = input.data(:,1);
time_raw = input.data(:,2);
time_raw_hrs = time_raw/3600;

figure(1)
plot(time_raw_hrs, OCP_raw)
xlabel('Time (hrs)')
ylabel('Open Circuit Potential (V)')
%{
idx_cut = find(time_raw_hrs>=48, 1);
OCP_cut = OCP_raw(idx_cut:end);
time_cut = time_raw_hrs(idx_cut:end);

figure(2)
plot(time_cut, OCP_cut)
xlabel('Time (hrs)')
ylabel('Open Circuit Potential (V)')
%}
%idx_cut = find(time_raw_hrs<=0.8, 1,'Last');
time_edit_hrs = time_raw_hrs;
OCP_edit = OCP_raw;
time_edit_hrs(OCP_edit<0.06)=[];
OCP_edit(OCP_edit<0.06)=[];
p_fit_all = polyfit(time_edit_hrs,OCP_edit,1);
disp(p_fit_all)

figure(3)
hold on
plot(time_edit_hrs, OCP_edit,'LineWidth',2)
plot(time_edit_hrs, polyval(p_fit_all,time_edit_hrs),'LineWidth',2)
hold off
legend('Raw Data',sprintf('Linear Fit - m = %0.2fmV/hr',p_fit_all(1)*1e3))
xlabel('Time (hrs)')
ylabel('Open Circuit Potential (V)')
set(gca,'FontSize',12)
set(gca,'LineWidth',1)

%{
idx_fit_1 = find(time_raw_hrs>=55.44, 1);
p_fit_1 = polyfit(time_raw_hrs(idx_fit:idx_fit_1),OCP_raw(idx_fit:idx_fit_1),1);
disp(p_fit_1)

idx_fit_2 = find(time_raw_hrs>=67.78, 1);
p_fit_2 = polyfit(time_raw_hrs(idx_fit_1:idx_fit_2),OCP_raw(idx_fit_1:idx_fit_2),1);
disp(p_fit_2)

p_fit_3 = polyfit(time_raw_hrs(idx_fit_2+1:end),OCP_raw(idx_fit_2+1:end),1);
disp(p_fit_3)
%}
%{
figure(4)
hold on
plot(time_cut, OCP_cut)
plot(time_raw_hrs(idx_fit:idx_fit_1), polyval(p_fit_1,time_raw_hrs(idx_fit:idx_fit_1)))
plot(time_raw_hrs(idx_fit_1:idx_fit_2), polyval(p_fit_2,time_raw_hrs(idx_fit_1:idx_fit_2)))
plot(time_raw_hrs(idx_fit_2+1:end), polyval(p_fit_3,time_raw_hrs(idx_fit_2+1:end)))
hold off
legend('Raw Data','Linear Fit - m = 2.7mV/hr', 'Linear Fit - m = -1.4mV/hr','Linear Fit - m = -0.4mV/hr')
xlabel('Time (hrs)')
ylabel('Open Circuit Potential (V)')

LineWidth = 1.5;
figure(5)
hold on
plot(time_cut, OCP_cut)
plot(time_cut, polyval(p_fit_all,time_cut),'LineWidth',LineWidth)
plot(time_raw_hrs(idx_fit:idx_fit_1), polyval(p_fit_1,time_raw_hrs(idx_fit:idx_fit_1)),'LineWidth',LineWidth)
plot(time_raw_hrs(idx_fit_1:idx_fit_2), polyval(p_fit_2,time_raw_hrs(idx_fit_1:idx_fit_2)),'LineWidth',LineWidth)
plot(time_raw_hrs(idx_fit_2+1:end), polyval(p_fit_3,time_raw_hrs(idx_fit_2+1:end)),'LineWidth',LineWidth)
hold off
legend('Raw Data','Linear Fit - m = 0.1mV/hr','Linear Fit - m = 2.7mV/hr', 'Linear Fit - m = -1.4mV/hr','Linear Fit - m = -0.4mV/hr')
xlabel('Time (hrs)')
ylabel('Open Circuit Potential (V)')

%% Porous RE
input2 = importdata('TAHU_OCP_500mMKCl.txt');
OCP_raw2 = input2.data(:,1);
time_raw2 = input2.data(:,2);
time_raw_hrs2 = time_raw2/3600;

figure(6)
plot(time_raw_hrs2, OCP_raw2)
xlabel('Time (hrs)')
ylabel('Open Circuit Potential (V)')

idx_fit2 = find(time_raw_hrs2>=0.81, 1);
p_fit_all2 = polyfit(time_raw_hrs2(1:idx_fit2),OCP_raw2(1:idx_fit2),1);
disp(p_fit_all2)

figure(7)
hold on
plot(time_raw_hrs2(1:end-10), OCP_raw2(1:end-10))
plot(time_raw_hrs2(1:idx_fit2), polyval(p_fit_all2,time_raw_hrs2(1:idx_fit2)))
hold off
legend('Raw Data','Linear Fit - m = -3.0mV/hr')
xlabel('Time (hrs)')
ylabel('Open Circuit Potential (V)')
%}
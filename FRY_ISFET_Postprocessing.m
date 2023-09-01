clear all

VgSweep = importdata('FRY12_ISFET1_VgSweep_pH7.xls');
Vg = VgSweep.data.Run321(:,9);
Vs = VgSweep.data.Run321(:,7);


figure(1)
clf
hold on
plot(Vg,Vs,'-b','LineWidth',2)
plot([0],[-1],'-r','LineWidth',2)
hold off
xlim([0 6])
ylim([-0.3 1.5])
legend('Experimental','Simulated')
xlabel('Gate Voltage (V)')
ylabel('Source Voltage (V)')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',14)


tSweep = importdata('FRY_Conv3.xls');
t = tSweep.data.Run197(:,1);
Vs_t = tSweep.data.Run197(:,7);
%}
%{
tSweep = importdata('FRY_Conv1.xls');
t = tSweep.data.Run193(:,1);
Vs_t = tSweep.data.Run193(:,7);
%}
%{
tSweep = importdata('FRY_Conv3.xls');
t = tSweep.data.Run195(:,1);
Vs_t = tSweep.data.Run195(:,7);
%}

%Conv3
a = 1.041;
b = 66.88;
c = 0.7786;
%}

%{
%Conv1
a=1.015;
b=97.06;
c=0.7958;
%}
%{
%Conv2
a=0.8407;
b=85.48;
c=0.5276;
%}
y_fit = a*(exp(-(t/b))-1)+c;

Vs_cor = Vs_t - y_fit;
h_cor = figure(3);
clf
plot(t(80:140),Vs_cor(80:140),'-g','LineWidth',2)
xlabel('Time (s)')
ylabel('Drift Corrected(V)')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',14)
shg

h_t = figure(2);
clf
hold on
plot(t,Vs_t)
plot(t,y_fit)
plot([0],[0],'-g','LineWidth',2)
hold off
xlabel('Time (s)')
ylabel('Source Voltage (V)')
legend('Raw data','Drift','Compensated','Location','southwest')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',14)
shg

[h_m h_i]=inset(h_t,h_cor);
set(gca,'LineWidth',1.5)
set(gca,'FontSize',14)
%legend('Data','Drift','Location','southwest')
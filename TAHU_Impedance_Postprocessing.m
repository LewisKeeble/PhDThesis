input = importdata('TAHU4_2_PotentiostaticACImpedance_5mVExcitation.csv');
f_agcl = input.data(:,1);
Z_agcl = input.data(:,4);
Z_prime_agcl = input.data(:,2);
Z_primeprime_agcl = input.data(:,3);
Phase_agcl = input.data(:,5);

plot_idx_agcl = find(f_agcl>=10);

figure(1)
plot(Z_prime_agcl(plot_idx_agcl),-Z_primeprime_agcl(plot_idx_agcl))
xlabel('Z''')
ylabel('-Z''''')

input = importdata('tahu4_2_enigce_potentiostaticacimpedance_5mvexcitation.csv');
f_enig = input.data(:,1);
Z_enig = input.data(:,4);
Z_prime_enig = input.data(:,2);
Z_primeprime_enig = input.data(:,3);
Phase_enig = input.data(:,5);

plot_idx_enig = find(f_enig>=10);

figure(2)
plot(Z_prime_enig(plot_idx_enig),-Z_primeprime_enig(plot_idx_enig))
xlabel('Z''')
ylabel('-Z''''')

figure(3)
hold on
plot(Z_prime_agcl(plot_idx_agcl),-Z_primeprime_agcl(plot_idx_agcl),'LineWidth',2)
plot(Z_prime_enig(plot_idx_enig),-Z_primeprime_enig(plot_idx_enig),'LineWidth',2)
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
hold off
xlabel('Z'' (\Omega)')
ylabel('-Z'''' (\Omega)')
legend('Ag/AgCl','ENIG')

figure(4)
hold on
plot(f_agcl(plot_idx_agcl),abs(Z_agcl(plot_idx_agcl)),'LineWidth',2)
plot(f_enig(plot_idx_enig),abs(Z_enig(plot_idx_enig)),'LineWidth',2)
set(gca,'XScale','log')
set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)
hold off
xlabel('Frequency (Hz)')
ylabel('Z (\Omega)')
legend('Ag/AgCl','ENIG')

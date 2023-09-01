in = importdata('Spectra_Au_1.txt');
E = in.data(:,1);
signal = in.data(:,2);

Estart = 0.75;
Estart_idx = find(E>=Estart,1,'first');

Eend = 3;
Eend_idx = find(E<=Eend,1,'last');

signal = movmean(log10(signal),5);

elements=['Si', 'Al', 'Zn', 'Ni', 'Au', 'Ag', 'Cl'];
lines_1 = [1.739, 1.486,1.012, 0.851, 2.12, 2.984, 2.621];
lines_2 = ['n\a','n\a',8.63, 7.471, 9.712, 'n\a', 'n\a'];

figure(1)
hold on
plot(E(Estart_idx:Eend_idx),signal(Estart_idx:Eend_idx),'LineWidth',4)
for i=[1,2,4,5]
    xline(lines_1(i),'--k','LineWidth',3)
end
hold off
xlim([Estart Eend])
ylim([min(signal(Estart_idx:Eend_idx)),max(signal(Estart_idx:Eend_idx))])
set(gca,'FontSize',18)
set(gca,'LineWidth',3)
xlabel('Energy (eV)')
ylabel('log_1_0(Counts)')

in = importdata('Spectra_Ag_1.txt');
E = in.data(:,1);
signal = in.data(:,2);

Estart = 0.75;
Estart_idx = find(E>=Estart,1,'first');

Eend = 4;
Eend_idx = find(E<=Eend,1,'last');

signal = movmean(log10(signal),5);

figure(2)
hold on
plot(E(Estart_idx:Eend_idx),signal(Estart_idx:Eend_idx),'LineWidth',4)
for i=[6,7]
    xline(lines_1(i),'--k','LineWidth',3)
end
hold off
xlim([Estart Eend])
ylim([min(signal(Estart_idx:Eend_idx)),max(signal(Estart_idx:Eend_idx))])
set(gca,'FontSize',18,'FontWeight','bold')
set(gca,'LineWidth',3)
xlabel('Energy (eV)')
ylabel('log_1_0(Counts)')


in = importdata('Spectra_Ag_2.txt');
E = in.data(:,1);
signal = in.data(:,2);

Estart = 0.75;
Estart_idx = find(E>=Estart,1,'first');

Eend = 4;
Eend_idx = find(E<=Eend,1,'last');

signal = movmean(log10(signal),5);

figure(3)
hold on
plot(E(Estart_idx:Eend_idx),signal(Estart_idx:Eend_idx),'LineWidth',4)
for i=[6,7]
    xline(lines_1(i),'--k','LineWidth',3)
end
hold off
xlim([Estart Eend])
ylim([min(signal(Estart_idx:Eend_idx)),max(signal(Estart_idx:Eend_idx))])
set(gca,'FontSize',18,'FontWeight','bold')
set(gca,'LineWidth',3)
xlabel('Energy (eV)')
ylabel('log_1_0(Counts)')
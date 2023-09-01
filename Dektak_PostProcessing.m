clear

Profile_ENIG = importdata('TAHU2_1_ENIG_Edited.csv');
Profile_Ag = importdata('TAHU2_1_WEAgPlated_Edited.csv');
Profile_Bare = importdata('TAHU2_1_Bare_Edited.csv');

idx_start_ENIG = 2855;
idx_start_Ag = 2726;
idx_start_Bare = 2624;

idx_end_ENIG = 48170;
idx_end_Ag = 48110;
idx_end_Bare = 58540;
bare_offset = 17.5;

figure(1)
hold on
plot(Profile_ENIG(idx_start_ENIG:idx_end_ENIG,1)-Profile_ENIG(idx_start_ENIG,1),Profile_ENIG(idx_start_ENIG:idx_end_ENIG,2)-Profile_ENIG(idx_end_ENIG,2))
plot(Profile_Ag(idx_start_Ag:idx_end_Ag,1)-Profile_Ag(idx_start_Ag,1),Profile_Ag(idx_start_Ag:idx_end_Ag,2)-Profile_Ag(idx_end_Ag,2))
plot(Profile_Bare(idx_start_Bare:idx_end_Bare,1)-Profile_Bare(idx_start_Bare,1)+bare_offset,Profile_Bare(idx_start_Bare:idx_end_Bare,2)-Profile_Bare(idx_end_Bare,2))
hold off
legend('ENIG','Ag','Bare')
shg

xlabel('Distance (\mum)')
ylabel('Height (\mum)')
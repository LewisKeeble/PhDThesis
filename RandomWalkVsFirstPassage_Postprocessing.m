
clear
reflecting_flag = 1;
ceiling_flag = 1;
if reflecting_flag == 1

    RW_1 = importdata('RandomWalkVsFirstPassage_BCs_RW_Ref_1.000000e-01z_NoCeiling');
    RW_2 = importdata('RandomWalkVsFirstPassage_BCs_RW_Ref_5.000000e-01z_NoCeiling');
    RW_3 = importdata('RandomWalkVsFirstPassage_BCs_RW_Ref_9.000000e-01z_NoCeiling');

    FP_1 = importdata('RandomWalkVsFirstPassage_BCs_FP_1.000000e-01z_2');
    FP_1 = FP_1(:,FP_1(4,:)<2400);
    FP_2 = importdata('RandomWalkVsFirstPassage_BCs_FP_5.000000e-01z_2');
    FP_2 = FP_2(:,FP_2(4,:)<2400);
    FP_3 = importdata('RandomWalkVsFirstPassage_BCs_FP_9.000000e-01z_2');
    FP_3 = FP_3(:,FP_3(4,:)<2400);

    if ceiling_flag == 1
        RC_ref_1 = importdata('RandomWalkVsFirstPassage//RandomWalkVsFirstPassage_BCs_RW_1.000000e-01z_2');
        RC_ref_2 = importdata('RandomWalkVsFirstPassage//RandomWalkVsFirstPassage_BCs_RW_5.000000e-01z_2');
        RC_ref_3 = importdata('RandomWalkVsFirstPassage//RandomWalkVsFirstPassage_BCs_RW_9.000000e-01z_2');
        RC_abs_1 = importdata('RWVsFP_RW_Ref_AbsCeiling_1.000000e-01z');
        RC_abs_2 = importdata('RWVsFP_RW_Ref_AbsCeiling_5.000000e-01z');
        RC_abs_3 = importdata('RWVsFP_RW_Ref_AbsCeiling_9.000000e-01z');

    end
else

    RW_1 = importdata('RandomWalkVsFirstPassage_BCs_RW_Abs_2_1.000000e-01z');
    RW_2 = importdata('RandomWalkVsFirstPassage_BCs_RW_Abs_2_5.000000e-01z');
    RW_3 = importdata('RandomWalkVsFirstPassage_BCs_RW_Abs_2_9.000000e-01z');

    if ceiling_flag == 1
        RC_abs_1 = importdata('RandomWalkVsFirstPassage_BCs_RW_Abs_NoCeiling_2_1.000000e-01z');
        RC_abs_2 = importdata('RandomWalkVsFirstPassage_BCs_RW_Abs_NoCeiling_2_5.000000e-01z');
        RC_abs_3 = importdata('RandomWalkVsFirstPassage_BCs_RW_Abs_NoCeiling_2_9.000000e-01z');
    end

    FP_1 = importdata('RandomWalkVsFirstPassage_BCs_FP_Abs_1.000000e-01z_2');
    FP_1 = FP_1(:,FP_1(4,:)<2400);
    FP_2 = importdata('RandomWalkVsFirstPassage_BCs_FP_Abs_5.000000e-01z_2');
    FP_2 = FP_2(:,FP_2(4,:)<2400);
    FP_3 = importdata('RandomWalkVsFirstPassage_BCs_FP_Abs_9.000000e-01z_2');
    FP_3 = FP_3(:,FP_3(4,:)<2400);

end
figure(1)
FP_h_1 = histfit(FP_1(1,:),100);

figure(2)
FP_h_2 = histfit(FP_2(1,:),100);

figure(3)
FP_h_3 = histfit(FP_3(1,:),100);

figure(4)
RW_h_1 = histfit(RW_1(1,:),100);

figure(5)
RW_h_2 = histfit(RW_2(1,:),100);

figure(6)
RW_h_3 = histfit(RW_3(1,:),100);


FP_dist_1 = fitdist(FP_h_1(2).XData','Normal');
FP_dist_2 = fitdist(FP_h_2(2).XData','Normal');
FP_dist_3 = fitdist(FP_h_3(2).XData','Normal');
RW_dist_1 = fitdist(RW_h_1(2).XData','Normal');
RW_dist_2 = fitdist(RW_h_2(2).XData','Normal');
RW_dist_3 = fitdist(RW_h_3(2).XData','Normal');

figure(7)
hold on
plot([0],'-k')
plot([0],'--k')
plot(FP_h_1(2).XData,FP_h_1(2).YData,'-r')
plot(FP_h_2(2).XData,FP_h_2(2).YData,'-b')
plot(FP_h_3(2).XData,FP_h_3(2).YData,'-g')
plot(RW_h_1(2).XData,RW_h_1(2).YData,'--r')
plot(RW_h_2(2).XData,RW_h_2(2).YData,'--b')
plot(RW_h_3(2).XData,RW_h_3(2).YData,'--g')
hold off
xlabel('Time (s)')
ylabel('x (m)')
legend('FP','RW')
xlim([-1 3]*1e-3);
ylim([0 2500]);

RW_mu = [RW_dist_1.mu,RW_dist_2.mu,RW_dist_3.mu];
FP_mu = [FP_dist_1.mu,FP_dist_2.mu,FP_dist_3.mu];

RW_sigma = [RW_dist_1.sigma,RW_dist_2.sigma,RW_dist_3.sigma];
FP_sigma = [FP_dist_1.sigma,FP_dist_2.sigma,FP_dist_3.sigma];

t = 0:1:2400;

FP_1_CS = zeros(1,length(t));
for i = t+1
FP_1_CS(i) = length(find(FP_1(4,:)<i));
end
FP_2_CS = zeros(1,length(t));
for i = t+1
FP_2_CS(i) = length(find(FP_2(4,:)<i));
end
FP_3_CS = zeros(1,length(t));
for i = t+1
FP_3_CS(i) = length(find(FP_3(4,:)<i));
end

RW_1_CS = zeros(1,length(t));
for i = t+1
RW_1_CS(i) = length(find(RW_1(4,:)<i));
end
RW_2_CS = zeros(1,length(t));
for i = t+1
RW_2_CS(i) = length(find(RW_2(4,:)<i));
end
RW_3_CS = zeros(1,length(t));
for i = t+1
RW_3_CS(i) = length(find(RW_3(4,:)<i));
end
if ceiling_flag == 1
    RC_abs_1_CS = zeros(1,length(t));
    for i = t+1
    RC_abs_1_CS(i) = length(find(RC_abs_1(4,:)<i));
    end
    RC_abs_2_CS = zeros(1,length(t));
    for i = t+1
    RC_abs_2_CS(i) = length(find(RC_abs_2(4,:)<i));
    end
    RC_abs_3_CS = zeros(1,length(t));
    for i = t+1
    RC_abs_3_CS(i) = length(find(RC_abs_3(4,:)<i));
    end
    RC_ref_1_CS = zeros(1,length(t));
    for i = t+1
    RC_ref_1_CS(i) = length(find(RC_ref_1(4,:)<i));
    end
    RC_ref_2_CS = zeros(1,length(t));
    for i = t+1
    RC_ref_2_CS(i) = length(find(RC_ref_2(4,:)<i));
    end
    RC_ref_3_CS = zeros(1,length(t));
    for i = t+1
    RC_ref_3_CS(i) = length(find(RC_ref_3(4,:)<i));
    end
end
figure(8)
hold on
plot([0],'-k')
plot([0],'--k')
plot([0],'-.k')
plot(FP_1_CS,'-','Color',[0.75 0 0])
plot(RW_1_CS,'--','Color',[1 0.25 0.25])
plot(FP_2_CS,'-','Color',[0 0 0.75])
plot(RW_2_CS,'--','Color',[0.25 0.25 1])
plot(FP_3_CS,'-','Color',[0 0.75 0])
plot(RW_3_CS,'--','Color',[0.25 1 0.25])
if ceiling_flag == 1
    plot(RC_ref_1_CS,'--r')
    plot(RC_ref_2_CS,'--b')
    plot(RC_ref_3_CS,'--g')
    plot(RC_abs_1_CS,'-.r')
    plot(RC_abs_2_CS,'-.b')
    plot(RC_abs_3_CS,'-.g')        
end
xlabel('Time (s)')
ylabel('Cumulative protons absorbed at z=0')
%legend('FP','RW - Reflecting Ceiling','RW - Absorbing Ceiling')
set(gca,'FontSize',12)
set(gca,'LineWidth',1.5)
hold off

shg
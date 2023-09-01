clear
t = 0:1:2400;

FP = load('RandomWalkVsFirstPassage_NoBCs_FP_4.mat');
FP = FP.protons_final;

FP_CS = zeros(1,length(t));
for i = 1:length(t)
FP_CS(i) = length(find(FP(4,:)<i));
end

RW = load('RandomWalkVsFirstPassage_NoBCs_RW_2.mat');
RW = RW.protons_final;

RW_CS = zeros(1,length(t));
for i = 1:length(t)
RW_CS(i) = length(find(RW(4,:)<i));
end


figure(1)
clf
hold on
plot(t,FP_CS)
plot(t,RW_CS)
hold off
shg

FP_x = FP(:,FP(1,:)>-1.1e-3&FP(1,:)<1.1e-3);
RW_x = RW(:,RW(1,:)>-1.1e-3&RW(1,:)<1.1e-3);

figure(2)
FP_h_1 = histfit(FP_x(1,:),100);

figure(3)
RW_h_1 = histfit(RW_x(1,:),100);

FP_dist_1 = fitdist(FP_h_1(2).XData','Normal');
RW_dist_1 = fitdist(RW_h_1(2).XData','Normal');

figure(7)
hold on
plot(FP_h_1(2).XData,FP_h_1(2).YData,'-r')
plot(RW_h_1(2).XData,RW_h_1(2).YData,'--b')
hold off
legend('FP','RW')
xlabel('Time (s)')
ylabel('Protons absorbed at z=0 plane')

RW_mu = RW_dist_1.mu;
FP_mu = FP_dist_1.mu;

RW_sigma = RW_dist_1.sigma;
FP_sigma = FP_dist_1.sigma;

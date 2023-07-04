t = 1:2401;
protons_released_FEM = squeeze(pH_FEM(:,2,:)*120e-6*120e-6*100e-6*6.02e23);
del_surfconc_FEM = (surfconc_FEM(:,2:4,:)-surfconc_FEM(1,2:4,:))*6.02e23*10e-6*10e-6;
sensor_nH_FEM = squeeze(del_surfconc_FEM(:,2,:)+del_surfconc_FEM(:,3,:));

figure(1)
hold on
plot(squeeze(mean(cumsum(sensor_nH_RW,2)))./squeeze(max(mean(cumsum(sensor_nH_RW,2))))')
plot(squeeze(mean(cumsum(sensor_nH_FP,2)))./squeeze(max(mean(cumsum(sensor_nH_FP,2))))')
plot(sensor_nH_FEM./max(sensor_nH_FEM));  
legend('QD','FP','FEM')
hold off
shg

figure(2)
hold on
plot(cumsum(protons_released_RW(:,1)))
plot(protons_released_FP(:,1))
plot(protons_released_FEM(:,1))
%set(gca, 'YScale', 'log')
legend('QD','FP','FEM')
hold off
shg

plot_range = 1:4;

% Set color vector
cStart_Red = uisetcolor([1 0.8 0.8],'Select a start color');
cEnd_Red = uisetcolor([1 0.0 0],'Select an end color');
c_Red = interp1([1;length(plot_range)],[cStart_Red;cEnd_Red],(1:length(plot_range))');

% Set color vector
cStart_Blue = uisetcolor([0.8 0.8 1],'Select a start color');
cEnd_Blue = uisetcolor([0 0.0 1],'Select an end color');
c_Blue = interp1([1;length(plot_range)],[cStart_Blue;cEnd_Blue],(1:length(plot_range))');

cStart_Green = uisetcolor([0.8 1 0.8],'Select a start color');
cEnd_Green = uisetcolor([0 1 0],'Select an end color');
c_Green = interp1([1;length(plot_range)],[cStart_Green;cEnd_Green],(1:length(plot_range))');


figure(3)
hold on
legend_labels = strings(3,1);
legend_labels(3) = 'First Passage';
legend_labels(2) = 'First Passage';
legend_labels(1) = 'Diffusion';

for i = 1:length(plot_range)
   if plot_range(i) == 4
    plot(t, squeeze(mean(cumsum(sensor_nH_RW,2)))./squeeze(max(mean(cumsum(sensor_nH_RW,2))))', '--','Color',c_Red(i,:),'HandleVisibility','off');
    p = plot(t, squeeze(mean(cumsum(sensor_nH_FP,2)))./squeeze(max(mean(cumsum(sensor_nH_FP,2))))', '-x','Color',c_Blue(i,:),'HandleVisibility','off'); 
    q = plot(t, sensor_nH_FEM./max(sensor_nH_FEM), '-o','Color',c_Green(i,:),'HandleVisibility','off'); 
    p.MarkerIndices = 1:5:length(t);
    q.MarkerIndices = 1:5:length(t);
   elseif i ~= length(plot_range)
    plot(t, squeeze(mean(cumsum(sensor_nH_RW,2)))./squeeze(max(mean(cumsum(sensor_nH_RW,2))))', '--','Color',c_Red(i,:), 'HandleVisibility','off');
    plot(t, squeeze(mean(cumsum(sensor_nH_FP,2)))./squeeze(max(mean(cumsum(sensor_nH_FP,2))))', '-','Color',c_Blue(i,:), 'HandleVisibility','off');      
    plot(t, sensor_nH_FEM./max(sensor_nH_FEM), '-','Color',c_Green(i,:), 'HandleVisibility','off');  
   else
    plot(t, squeeze(mean(cumsum(sensor_nH_RW,2)))./squeeze(max(mean(cumsum(sensor_nH_RW,2))))', '--','Color',c_Red(i,:));
    plot(t, squeeze(mean(cumsum(sensor_nH_FP,2)))./squeeze(max(mean(cumsum(sensor_nH_FP,2))))', '-','Color',c_Blue(i,:)); 
    plot(t, sensor_nH_FEM./max(sensor_nH_FEM), '-','Color',c_Green(i,:));   
   end
   xlabel('t (s)')
   ylabel('Normalised mean array signal (# protons)')
   title('Benchmark signal') 
end
legend(legend_labels);
%###########################################################
% s_plot_expt_diagnostics.m
%###########################################################
% description:
%-----------------------------------------------------------
% This script plots experimental data of diagnostics.
%###########################################################


%author: Markus Markl
%created: 09.06.2021

PRINT = true;

addpath('~/Dokumente/plasma/code/libneo_old/matlab/Utility/xxaxis');

%shot = 33133;
%time = [2100, 2200, 2400, 2600, 2800, 3000];
shot = 34548;
%time =[1400, 1700, 1850, 2000, 2200, 2325, 2670, 2900];
[time, type] = get_shot_times(shot, '/temp/markl_m/ELMsuppression_in_hydrogen/DATA/INDEX/');
runname = [num2str(shot)];%, '_', num2str(time)];
plotname =  [runname, '_exptdiag'];

datapath = ['/temp/markl_m/ELMsuppression_in_hydrogen/DATA/ELMDIAG/', num2str(shot),'/'];
outpath = ['/temp/markl_m/'];

Ipolsoli = [num2str(shot), '_MAC_Ipolsoli.dat'];
IBl6 = [num2str(shot), '_MAW_IBl6.dat'];
ELMi = [num2str(shot), '_POT_ELMi-Han.dat'];

Ipolsoli_data = importdata([datapath, Ipolsoli]);
IBl6_data = importdata([datapath, IBl6]);
ELMi_data = importdata([datapath, ELMi]);

I_ymax = 2.1e4;
Da_ymax = 0.2;

fig = figure;

subplot(3,1,1)
plot(Ipolsoli_data(:,1), Ipolsoli_data(:,2), 'Color', 'blue');
title(['Shot ', num2str(shot), ': Divertor Current (MAC, Ipolsoli)']);
ylabel('I/A')
%xlabel('t/s')
%yl = ylim;
%ylim([0, yl(2)]);
ylim([0,I_ymax]);
yl = ylim;
xl = xlim;
xlim([0, xl(2)]);
hold on
for i= 1:numel(time)
	timeline = xline(time(i)/1000, '--');%, {[num2str(time(i))]});%, 'LineWidth', 2,);
	timeline.LineWidth = 2;
end
hold off

subplot(3,1,2)
plot(ELMi_data(:,1), ELMi_data(:,2), 'Color', 'red');
title(['Shot ', num2str(shot), ': D\alpha Radiation (POT, ELMi-Han)']);
ylabel('D\alpha /1')
%xlabel('t/s')
xlim(xl)
xlim([0, xl(2)]);
ylim([0, Da_ymax]);
hold on
for i= 1:numel(time)
	timeline = xline(time(i)/1000, '--');%, {[num2str(time(i))]});%, 'LineWidth', 2,);
	timeline.LineWidth = 2;
end
hold off

subplot(3,1,3)
plot(IBl6_data(:,1), IBl6_data(:,2), 'Color', 'k');
title(['Shot ', num2str(shot), ': Bl6 Coil Current (MAW, IBl6)']);
ylabel('I/kA')
xlabel('t/s')
xlim([0, xl(2)]);
hold on
for i= 1:numel(time)
	timeline = xline(time(i)/1000, '--');%, {[num2str(time(i))]});%, 'LineWidth', 2,);
	timeline.LineWidth = 2;
end
hold off

suptitle(' ');

if PRINT == true
	system(['mkdir -p ', outpath, 'exptdiag/']);
	print([outpath, 'exptdiag/', plotname, '.png'], '-dpng', '-r200');
	print([outpath, 'exptdiag/', plotname, '.svg'], '-dsvg');
end

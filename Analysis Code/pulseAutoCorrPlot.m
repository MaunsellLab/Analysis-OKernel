
% AutoCorrelation of White Noise Opto Pulses
staShift = 12; % offset of the sta peak. 
pulse = zeros(1,800);
pulse(376-staShift:400-staShift) = 1;
r = xcorr(pulse, pulse);

% plot
figure;
plot(-1*r(401+staShift:900+staShift), 'Color','k', 'LineWidth',2);
hold on; title('Pulse Autocorrelation');
plot([400 400], [-35 5], 'Color', 'k', 'LineStyle','--', 'LineWidth', 0.5);
ylabel('Normalized Power'); xlabel('Time Relative to Spike Onset (ms)');
set(gca, 'TickDir', 'out'); set(gca,'FontSize', 16);
box off; xticks([0 100 200 300 400 500]);
xticklabels({'-400','-300','-200','-100', '0', '100'})
yticks([-25 0 5]); ylim([-35 5]);
hold off;
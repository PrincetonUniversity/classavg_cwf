SNR=1/80;
fname=sprintf('/scratch/ARCHIVE_from_sdl6/tbhamre/classavg_paper_data/mah_hist_snr1by%d.mat',1/SNR); 
load(fname)
figure(1);
clf;
[counts1, binCenters1] = hist(ang_dm, 100);
[counts2, binCenters2] = hist(ang_df, 100);
plot(binCenters1, counts1/sum(counts1), 'b-', 'LineWidth',2);
hold on;
plot(binCenters2, counts2/sum(counts2), 'g-', 'LineWidth',2);
grid on;
% Put up legend.
legend1 = sprintf('Improved, Mean = %.3f', mean(ang_dm));
legend2 = sprintf('Initial, Mean = %.3f', mean(ang_df));
legend({legend1, legend2 }, 'Box', 'off','FontSize',18);
xlabel('Angular distance in degrees','FontSize',18)
ylabel('Probability Density Function','FontSize',18)

grid off
titlstr=sprintf('SNR=1/%d',1/SNR)
title(titlstr, 'FontSize',18)
fname=sprintf('fighist_snr1by%d.png',1/SNR)
fpath = '~/cwf_classavg/paper/';
print('-dpng',fullfile(fpath, fname)); % Save as vector graphics

fpath = sprintf('/scratch/ARCHIVE_from_sdl6/tbhamre/classavg_paper_data/mah_hist_snr1by%d.mat',1/SNR);
save(fpath,'ang_dm', 'ang_df')

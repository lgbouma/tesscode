noisedat = fitsread('noises.fits');

%sedat = fitsread('superearths_13x24.fits');
%etdat = fitsread('earths_13x24.fits');
%mndat = fitsread('mininep_13x24.fits');
%hzdat = fitsread('small_hz_13x24.fits');

%fitsout = [[imag], [npix], [noise], [shot_noise], [bknd_noise], [read_noise], [sys_noise], [satn], [diln]]

imag = noisedat(1,:);
npix = noisedat(2,:);
dil = noisedat(9,:);
tot = noisedat(3,:).*(1+dil);
shot = noisedat(4,:).*(1+dil);
bknd = noisedat(5,:).*(1+dil);
read = noisedat(6,:).*(1+dil);
sys = noisedat(7,:).*(1+dil);
sat = noisedat(8,:);


logsig = -5:0.1:-2.1;

magbins = (0:50)/3 + 5;

figure;
subplot(2,1,1)
% imagesc(magbins, logsig, etdat, [0 20]);
% colormap 'gray'
% colormap(flipud(colormap));
% set(gca, 'YDir', 'normal');
% c = colorbar;
% ylabel(c, 'Earths')
% hold on
plot(imag, log10(shot), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
hold on
plot(imag, log10(bknd), '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
plot(imag, log10(read), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
plot(imag, log10(sys), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
plot(imag, log10(tot), 'k-', 'LineWidth', 1)
xlabel('I_C', 'FontSize', 12)
ylabel('log_{10}(\sigma) [hr^{-1/2}]', 'FontSize', 12)
legend('Shot Noise', 'Sky Noise', 'Read Noise', 'Sys. Noise', 'Total');
%et_it = zeros(1,length(imag));
%for ii=1:length(imag)
%    et_it(ii) = sum(etdat(logsig>log10(noisedat(3,ii)),ii));
%end
%title([num2str(round(sum(et_it))) ' Earths']);
axis([5 18 -5 -1]);
set(gca, 'FontSize', 12);

subplot(2,1,2);
plot(imag, npix, 'k-', 'LineWidth', 1);
axis([5 18 0 20]);
xlabel('I_C', 'FontSize', 12)
ylabel('Pixels in Optimal Aperture', 'FontSize', 12);
set(gca, 'FontSize', 12);

% subplot(2,2,2);
% imagesc(magbins, logsig, sedat, [0 20]);
% colormap 'gray'
% colormap(flipud(colormap));
% set(gca, 'YDir', 'normal');
% c = colorbar;
% ylabel(c, 'Super-Earths')
% hold on
% plot(imag, log10(noisedat(4,:)), 'g-')
% plot(imag, log10(noisedat(5,:)), 'g-.')
% plot(imag, log10(noisedat(6,:)), 'g--')
% plot(imag, log10(noisedat(7,:)), 'g:')
% plot(imag, log10(noisedat(3,:)), 'r-', 'LineWidth', 1)
% xlabel('I_C')
% ylabel('log_{10}(\sigma) [hr^{-1/2}]')
% %legend('Shot Noise', 'Sky Noise', 'Read Noise', 'Sys. Noise', 'Total');
% se_it = zeros(1,length(imag));
% for ii=1:length(imag)
%     se_it(ii) = sum(sedat(logsig>log10(noisedat(3,ii)),ii));
% end
% title([num2str(round(sum(se_it))) ' Super-Earths']);
% axis([5 20 -5 -2]);
% 
% subplot(2,2,3);
% imagesc(magbins, logsig, mndat, [0 20]);
% colormap 'gray'
% colormap(flipud(colormap));
% set(gca, 'YDir', 'normal');
% c = colorbar;
% ylabel(c, 'Mini-Neptunes')
% hold on
% plot(imag, log10(noisedat(4,:)), 'g-')
% plot(imag, log10(noisedat(5,:)), 'g-.')
% plot(imag, log10(noisedat(6,:)), 'g--')
% plot(imag, log10(noisedat(7,:)), 'g:')
% plot(imag, log10(noisedat(3,:)), 'r-', 'LineWidth', 1)
% xlabel('I_C')
% ylabel('log_{10}(\sigma) [hr^{-1/2}]')
% %legend('Shot Noise', 'Sky Noise', 'Read Noise', 'Sys. Noise', 'Total');
% mn_it = zeros(1,length(imag));
% for ii=1:length(imag)
%     mn_it(ii) = sum(mndat(logsig>log10(noisedat(3,ii)),ii));
% end
% title([num2str(round(sum(mn_it))) ' Mini-Neptunes']);
% axis([5 20 -5 -2]);
% 
% subplot(2,2,4);
% imagesc(magbins, logsig, hzdat, [0 20]);
% colormap 'gray'
% colormap(flipud(colormap));
% set(gca, 'YDir', 'normal');
% c = colorbar;
% ylabel(c, 'Small HZ Planets')
% hold on
% plot(imag, log10(noisedat(4,:)), 'g-')
% plot(imag, log10(noisedat(5,:)), 'g-.')
% plot(imag, log10(noisedat(6,:)), 'g--')
% plot(imag, log10(noisedat(7,:)), 'g:')
% plot(imag, log10(noisedat(3,:)), 'r-', 'LineWidth', 1)
% xlabel('I_C')
% ylabel('log_{10}(\sigma) [hr^{-1/2}]')
% %legend('Shot Noise', 'Sky Noise', 'Read Noise', 'Sys. Noise', 'Total');
% hz_it = zeros(1,length(imag));
% for ii=1:length(imag)
%     hz_it(ii) = sum(hzdat(logsig>log10(noisedat(3,ii)),ii));
% end
% title([num2str(round(sum(hz_it))) ' Small HZ Planets']);
% axis([5 20 -5 -2]);


    
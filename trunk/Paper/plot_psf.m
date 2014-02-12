figure;
subplot(1,2,2);

dat1 = fitsread('../Background/besancon/splines.fits');
xplot = dat1(:,1:4);
yplot = dat1(:,5:8);

plot(xplot, yplot, 'LineWidth', 1);
set(gca, 'FontSize', 12);
legend('0\circ', '6\circ', '12\circ', '17\circ');
axis([0 7 0 11]);

xlabel('Distance [pixels]', 'Fontsize', 12);
ylabel('Rejection [mag]', 'Fontsize', 12);

subplot(1,2,1);
dat2 = fitsread('../ExpTimeCalc/frac24_105_4700k.fits');
xplot = 1:256;
yplot = squeeze(mean(mean(dat2,1),2));

semilogx(xplot, yplot, 'LineWidth', 1);
set(gca, 'FontSize', 12);
axis([1 100 0 1]);
xlabel('Pixels in Aperture', 'Fontsize', 12);
ylabel('Enclosed Flux Fraction', 'Fontsize', 12);
kopp_hz_in  = load('kopp_hz_in.txt');
kopp_hz_out = load('kopp_hz_out.txt');

dat = fitsread('../Report/mout_6_13x24_69.fits');
trial = dat(1,:);
rpla  = dat(11,:);

trial1 = find((trial==1)&(rpla<2));
teff  = dat(10,trial1);
spla  = dat(14,trial1);

teffplot = 3000:100:6000;

%kopp_sin_plot = interp1(kopp_hz_in(:,2), kopp_hz_in(:,1), teffplot, 'nearest');
%kopp_sout_plot = interp1(kopp_hz_out(:,2), kopp_hz_out(:,1), teffplot, 'nearest');

p = polyfit(kopp_hz_in(:,2)/1000, kopp_hz_in(:,1), 5);
kopp_sin_plot = polyval(p, teffplot/1000);
p = polyfit(kopp_hz_out(:,2)/1000, kopp_hz_out(:,1), 5);
kopp_sout_plot = polyval(p, teffplot/1000);

simp_sin_plot = 4.0*ones(size(teffplot));
simp_sout_plot = 0.25*ones(size(teffplot));

figure;plot(spla, teff, 'k+', 'LineWidth', 1);
hold on
plot(simp_sout_plot, teffplot, 'b:', 'LineWidth', 2);
plot(kopp_sout_plot, teffplot, 'b--', 'LineWidth', 1);
legend('TESS Detections', 'Petigura (2014) HZ', 'Kopparapu (2013) HZ');
%plot(kopp_hz_in(:,1), kopp_hz_in(:,2), 'r:');
%plot(kopp_hz_out(:,1), kopp_hz_out(:,2), 'b:');
plot(simp_sin_plot, teffplot, 'r:', 'LineWidth', 2);
plot(kopp_sin_plot, teffplot, 'r--', 'LineWidth', 1);
axis([0 6 3000 5000]);
set(gca, 'XDir', 'Reverse', 'FontSize', 12);
xlabel('S/S_0', 'FontSize', 12);
ylabel('T_{eff} [K]');
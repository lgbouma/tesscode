dat = fitsread('pointings.fits');
len = max(size(dat));
%plotind = ceil(rand(1,100000)*len);
plotind=1:len;
elon = dat(1,plotind);
elat = dat(2,plotind);
npnt = dat(3,plotind);
x = cosd(elon).*(90-abs(elat));
y = sind(elon).*(90-abs(elat));
%[xm ym] = meshgrid(x,y);
xplot = -90:0.2:90;
yplot = -90:0.2:90;
[X Y] = meshgrid(xplot, yplot);
%zplot = interp2(xm,ym,zm,xplot,yplot, 'linear', 0);
%scatter(x(plotind), y(plotind), 20, npnt(plotind), 'filled')
Z = griddata(x,y,npnt, X,Y);
%[xm ym zm] = meshgrid(x, y, npnt);
imagesc(xplot, yplot, Z);
axis([-95 95 -95 95]);
axis square
labels = [-90 -45 0 45 90];
set(gca, 'YDir', 'normal', 'Fontsize', 12, ...
 'XTick', labels, 'XTickLabel', labels, ...
 'YTick', labels, 'YTickLabel', labels);
colormap 'gray'
colormap(flipud(colormap));
c = colorbar;
set(c, 'Fontsize', 12);
ylabel(c, 'Number of Pointings', 'Fontsize', 12);
title('Ecliptic Polar Projection [Degrees]', 'Fontsize', 12);
hold on
theta = 0:0.001:2*pi;
polar(theta, 90*ones(size(theta)), 'k--');
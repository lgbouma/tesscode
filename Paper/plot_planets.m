dat = fitsread('planets.fits');

nstars = 1e6;

r = dat(1,:);
p = dat(2,:);

p_bound = [0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 418.0];
r_bound = [0.8, 1.25, 2.0, 4.0, 6.0, 22.0];

logp_bound = log10([0.8, 2.0, 3.4, 5.9, 10.0, 17.0, 29.0, 50.0, 85.0, 145.0, 418.0]);
logr_bound = log10([0.8, 1.25, 2.0, 4.0, 6.0, 22.0]);

p_label = [0.8, 2.0, 5.9, 17.0, 50.0, 145.0, 418.0];
r_label = [0.8, 1.25, 2.0, 4.0, 6.0, 22.0];

logp_label = log10([0.8, 2.0, 5.9, 17.0,  50.0, 145.0, 418.0]);
logr_label = log10([0.8, 1.25, 2.0, 4.0, 6.0, 22.0]);

logp = log10(p_bound);
logr = log10(r_bound);

dp = median(diff(logp))/6;
dr = median(diff(logr))/6;

logrgrid = (logr(1):dr:logr(end));
logpgrid = (logp(1):dp:logp(end));

rgrid = 10.^logrgrid;
pgrid = 10.^logpgrid;

ptick = interp1(logpgrid, 1:length(pgrid), logp_label, 'nearest', length(pgrid));
rtick = interp1(logrgrid, 1:length(rgrid), logr_label, 'nearest', length(rgrid));
% 
% hist2d = zeros(length(rgrid)-1, length(pgrid)-1);
% 
% for ii=1:(length(rgrid)-1)
% 	for jj=1:(length(pgrid)-1)
%         hist2d(ii,jj) = sum((r>rgrid(ii)) & (r<rgrid(ii+1)) & ...
%             (p>pgrid(jj)) & (p<pgrid(jj+1)));
%     end
% end
%      
% histr = sum(hist2d,2)/nstars/dr;
% histp = sum(hist2d,1)/nstars/dp;

figure;
subplot(2,2,1);
%imagesc(logpgrid(1:end), logrgrid(1:end), log10(hist2d));
imagesc(log10(hist2d));
ylabel('Radius [R_{\oplus}]', 'Fontsize', 12);
xlabel('Period [days]', 'Fontsize', 12);
title('Fressin Planet Distribution', 'FontSize', 12);
%set(gca, 'XTick', ptick-0.5*dp, 'XTickLabel', p_bound);
%set(gca, 'YTick', rtick-0.5*dr, 'YTickLabel', r_bound);
set(gca, 'XTick', ptick-0.5, 'XTickLabel', p_label);
set(gca, 'YTick', rtick-0.5, 'YTickLabel', r_label);
set(gca, 'YDir', 'normal', 'Fontsize', 12);
colormap 'gray';
colormap(flipud(colormap));

subplot(2,2,2);
b = barh(histr, 'histc');
%barh(histp,  'histc');
set(gca, 'YTick', rtick, 'YTickLabel', r_label, 'FontSize', 12);
axis([ 0 1.1*max(histr) 1 length(rgrid)]);
ylabel('Radius [R_{\oplus}]', 'Fontsize', 12);
xlabel('dN/dlog(R)');
title('Planet Occurrence', 'Fontsize', 12);

subplot(2,2,3);
bar(histp, 'histc');
set(gca, 'XTick', ptick, 'XTickLabel', p_label, 'FontSize', 12);
xlabel('Period [days]', 'Fontsize', 12);
ylabel('dN/dlog(P)');
title('Planet Occurrence', 'Fontsize', 12);
axis([1 length(pgrid) 0 1.1*max(histp)]);

dat = fitsread('planetonly.fits');

makecat = true;

trial = dat(1,:); %[eclip[det].trial], 
vmag  = dat(2,:); %[targets[detid].mag.v], 
imag  = dat(3,:); %[targets[detid].mag.ic], 
tmag  = dat(4,:); %[targets[detid].mag.t], 
jmag  = dat(5,:); %[targets[detid].mag.j], 
hmag  = dat(6,:); %[targets[detid].mag.j], 
kmag  = dat(7,:); %[targets[detid].mag.k], 
teff  = dat(8,:); %[targets[detid].teff], 
elon  = dat(9,:); %[eclip[det].coord.elon], 
elat  = dat(10,:); %[eclip[det].coord.elat], 
ra  = dat(11,:); %[eclip[det].coord.glon], 
dec  = dat(12,:); %[eclip[det].coord.glat], 
p     = dat(13,:); %[eclip[det].p], 
a     = dat(14,:); %[eclip[det].a], 
s     = dat(15,:); %[eclip[det].s], 
b     = dat(16,:); %[eclip[det].b], 
teff2 = dat(17,:); %[eclip[det].teff2], 
m2    = dat(18,:); %[eclip[det].m2], 
r2    = dat(19,:)/0.00917; %[eclip[det].r2], 
dep1  = dat(20,:); %[eclip[det].dep1], 
dur1  = dat(21,:); %[eclip[det].dur1], 
necl1 = dat(22,:); %[eclip[det].neclip_obs1], 
teff1 = dat(23,:); %[eclip[det].teff1], 
m1    = dat(24,:); %[eclip[det].m1], 
r1    = dat(25,:); %[eclip[det].r1], 
dep2  = dat(26,:); %[eclip[det].dep2], 
dur2  = dat(27,:); %[eclip[det].dur2], 
necl2 = dat(28,:); %[eclip[det].neclip_obs2], 
snr1  = dat(29,:); %[eclip[det].snreclp1], 
snrg1  = dat(30,:); %[eclip[det].snrgress1], 
snr2  = dat(31,:); %[eclip[det].snreclp2], 
snrg2  = dat(32,:); %[eclip[det].snrgress2], 
rvk   = dat(33,:); %[eclip[det].k], 
snrh  = dat(34,:); %[eclip[det].snrhr], 
starph = dat(35,:); %[eclip[det].star_ph],
bkph  = dat(36,:); %[eclip[det].bk_ph], 
zodi  = dat(37,:); %[eclip[det].zodi_ph], 
npix  = dat(38,:); %[eclip[det].npix],
dil   = dat(39,:); %[eclip[det].dil], 
ffi   = dat(40,:); %[targets[detid].ffi], 
npts  = dat(41,:); %[eclip[det].npointings] ,
sat   = dat(42,:); %[eclip[det].sat], 
fovr  = dat(43,:); %[eclip[det].coord.fov_r], 
eclass = dat(44,:);
hsep = dat(45,:);
icsys = dat(46,:);
tsys = dat(47,:);
jsys = dat(48,:);
censhift1 = dat(49,:);
censhift2 = dat(50,:);
cenerr1 = dat(51,:);
cenerr2 = dat(52,:);
bin   = dat(53,:); %[bins], 
binsep  = dat(54,:); %[targets[detid].companion.sep], 
bint  = dat(55,:); %[targets[targets[detid].companion.ind].mag.t]]


ntrials = max(trial);

earths = zeros(ntrials,1);
earthps = zeros(ntrials,1);
earthffi = zeros(ntrials,1);
speths = zeros(ntrials,1);
spethps = zeros(ntrials,1);
spethffi = zeros(ntrials,1);
smneps = zeros(ntrials,1);
smnepps = zeros(ntrials,1);
smnepffi = zeros(ntrials,1);
smplas = zeros(ntrials,1);
hzplas = zeros(ntrials,1);
hzplasps = zeros(ntrials,1);
hzeclp = zeros(ntrials,1);
hzeclpps = zeros(ntrials,1);
giants = zeros(ntrials,1);
giantps = zeros(ntrials,1);
coolgiants = zeros(ntrials,1);
binplans = zeros(ntrials,1);
brights = zeros(ntrials,1);
aroundm = zeros(ntrials,1);

for ii=1:ntrials
    earths(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 1.25));
    aroundm(ii) = sum((eclass==1) & (trial==ii) & (teff<3500) & (jmag<10));
    brights(ii) = sum((eclass==1) & (trial==ii)&(r2>2)&(r2<4)&(imag<12)& ~ffi);
    brights(ii) = sum((eclass==1) & (trial==ii)&(r2<4)&(imag<12)&(rvk>1)& ~ffi);
    earthps(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 1.25) & ~ffi);
    earthffi(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 1.25) & ffi);
    speths(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0) & (r2 > 1.25));
    spethps(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0) & (r2 > 1.25) & ~ffi);
    spethffi(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0) & (r2 > 1.25) & ffi);
    smneps(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 4) & (r2 > 2.0));
    smnepps(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 4) & (r2 > 2.0) & ~ffi);
    smnepffi(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 4) & (r2 > 2.0) & ffi);
    smplas(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0));
    giants(ii) = sum((eclass==1) & (trial==ii) & (r2 > 4.0));
    giantps(ii) = sum((eclass==1) & (trial==ii) & (r2 > 4.0) & ~ffi);
    coolgiants(ii) = sum((eclass==1) & (trial==ii) & (r2 > 4.0) & (s<150) & (vmag<12) & ffi);
    hzplas(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0) & ...
        (s <= 2) & (s >= 0.5));
    hzplasps(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0) & ...
        (s <= 2) & (s >= 0.5) & ~ffi);
    hzeclp(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0) & ...
        (s <= 2) & (s >= 0.5) & ...
        (abs(elat) >= 75));
    hzeclpps(ii) = sum((eclass==1) & (trial==ii) & (r2 <= 2.0) & ...
        (s <= 2) & (s >= 0.5) & ...
        (abs(elat) >= 75) & ~ffi);
    binplans(ii) = sum((eclass==1) & (trial==ii)&(bin>0));
end

[   num2str(mean(earths),3) ' $\pm$ ' num2str(std(earths),2) ' & '...
    num2str(mean(speths),3) ' $\pm$ ' num2str(std(speths),2) ' & '...
    num2str(mean(smplas),3) ' $\pm$ ' num2str(std(smplas),2) ' & '...
    num2str(mean(smneps),3) ' $\pm$ ' num2str(std(smneps),2) ' & '...
    num2str(mean(hzplas),3) ' $\pm$ ' num2str(std(hzplas),2) ' & '...
    num2str(mean(hzeclp),3) ' $\pm$ ' num2str(std(hzeclp),2) ' \\ ']

plas = find(eclass==1);
[y, ind] = sort(r2(plas));
rsort = plas(ind);
rsort = rsort(1:20);

if makecat
   fid = fopen('test.cat', 'w');
   outdat =  [ra(rsort)' dec(rsort)' r2(rsort)' p(rsort)' s(rsort)' teff1(rsort)' vmag(rsort)' imag(rsort)' jmag(rsort)' kmag(rsort)']';
   fprintf(fid, '%6.3f & %6.3f & %3.2f & %3.2f & %4.2f & %4.0f & %4.2f & %4.2f & %4.2f & %4.2f \\\\ \n', outdat);
   fclose(fid);
   %dlmwrite('test.cat', [ra(rsort)' dec(rsort)' r2(rsort)' p(rsort)' s(rsort)' teff1(rsort)' vmag(rsort)' imag(rsort)' jmag(rsort)' kmag(rsort)'], 'delimiter', ' & ', 'precision', 4);    
end


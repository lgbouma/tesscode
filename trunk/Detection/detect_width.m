% 1-minute cadence to data, so all times in minutes
% Campaign duration:
D = round(27*24*60);

% Transit width for circular orbit:
%W = (1.3*60)*(P/(24*60)).^(1/3);
Wmin = floor(78*(4.6/24)^(1/3));
Wmax = ceil(78*(27/2)^(1/3));
W = Wmin:Wmax;

% Number of transit widths:
NW = floor(D./W);

% Period in minutes
P = 24*60*(W/78).^3;

% Period in width range
Pwmin = floor(P./W);
Pwmax = ceil(P./W);

% Number of guaranteed transits
NTmin = ceil(NW./Pwmin);
NTmax = ceil(NW./Pwmax);

%Number of stars
Nstar = 5e2;
% Mean of transit trials:
TSM = zeros(1,sum(NW));
tsmind = cumsum(NW);
% Transit sub-sampling
NSUB = 2;
% Save threshold
NSIG = 4;
% Number of averaged points
%L = W.*NT;
% STD in phase-folded light curve
%sig_sum = 1./sqrt(L);
% Events out
hisigs = [];

for jj=1:Nstar
  % Normally-distributed timeseries:
  TS = randn(D,1);
  for kk=1:NSUB
    for ii=1:length(W)  
      % circ-shift for sub-sampling
      if kk>1
        os = round((kk-1)*W(ii)/NSUB);
        % Truncate the lightcurve to a multiple of the widths
        TST = [TS((os+1):(W(ii)*NW(ii))); TS(1:os)];
      else
        TST = TS(1:(W(ii)*NW(ii)));
      end
      % Width-fold the lightcurve
      % Take mean over each possible transit
      TSW = mean(reshape(TST,NW(ii),W(ii)),2);
      cts = ones(size(TSW));
      
      % Zero-pad
      TSWmin = [TSW; zeros(NTmin(ii)*Pwmin(ii)-length(TSW),1)];
      ctsmin = [cts; zeros(NTmin(ii)*Pwmin(ii)-length(TSW),1)];
      TSWmax = [TSW; zeros(NTmax(ii)*Pwmax(ii)-length(TSW),1)];
      ctsmax = [cts; zeros(NTmax(ii)*Pwmax(ii)-length(TSW),1)];
      csmin = sum(reshape(ctsmin, Pwmin(ii), NTmin(ii)),2);
      csmax = sum(reshape(ctsmax, Pwmax(ii), NTmax(ii)),2);
      
      % Modified mean
      TSPmin = sum(reshape(TSWmin, Pwmin(ii), NTmin(ii)),2)./csmin;          
      TSPmax = sum(reshape(TSWmax, Pwmax(ii), NTmax(ii)),2)./csmax;
      
      % Take the ntransit>1 cases 
      TSMmin = TSPmin(csmin>1).*sqrt(W(ii)*csmin(csmin>1));
      TSMmax = TSPmax(csmax>1).*sqrt(W(ii)*csmax(csmax>1));
      minevent = (TSPmin>NSIG);
      maxevent = (TSPmax>NSIG);
      nminevents = sum(minevent);
      nmaxevents = sum(maxevent);
      %display([num2str(nevents) ' events on trial ' num2str(jj)]);
      if (nminevents>0)
        hisigs = [hisigs TSMmin(minevent)];
        %display([num2str(nminevents) ' events on trial ' num2str(jj)]);
      end 
      if (nmaxevents>0)
        hisigs = [hisigs TSMmax(maxevent)];
      end 
    end
  end
end

save 'hisigs.mat' hisigs
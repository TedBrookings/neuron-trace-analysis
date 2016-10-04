function bursts = GetBurstQuantification( burstStartTimes, burstStopTimes, ...
                                          spikes, t )
  if ~isrow( burstStartTimes ), burstStartTimes = burstStartTimes'; end
  if ~isrow( burstStopTimes ), burstStopTimes = burstStopTimes'; end
  numBursts = numel( burstStartTimes );
  burstSpikeInds = cell( 1, numBursts );
  nonBurstSpikeInds = cell( 1, numBursts + 1 );

  t2 = -Inf;
  for n = 1:numBursts
    t1 = t2; t2 = burstStartTimes(n);
    % note <, not <=
    nonBurstSpikeInds{n} = find( t1 < spikes.maxV.t & spikes.maxV.t < t2 );
    t1 = t2; t2 = burstStopTimes(n);
    % note <=, not <
    burstSpikeInds{n} = find( t1 <= spikes.maxV.t & spikes.maxV.t <= t2 );
  end
  t1 = t2; t2 = t(end);
  nonBurstSpikeInds{numBursts+1} ...
    = find( t1 <= spikes.maxV.t & spikes.maxV.t <= t2 );

  burstStartTimes = burstStartTimes .* 1e-3;
  burstStopTimes = burstStopTimes .* 1e-3;
  burstDurations = burstStopTimes - burstStartTimes;
  interBurstIntervals = burstStopTimes(2:end) - burstStartTimes(1:end-1);
  if isempty( burstStartTimes )
    allNonBurstIntervals = (t(end) - t(1)) * 1e-3;
  else
    lastInterval = t(end) * 1e-3 - burstStopTimes(end);
    allNonBurstIntervals = [ burstStartTimes(1), interBurstIntervals, ...
                             lastInterval ];
  end
  burstPeriods = diff( burstStartTimes );
  burstDutyCycle = burstDurations(1:end-1) ./ burstPeriods;
  burstRates = 1.0 ./ burstPeriods;
  numSpikesPerBurst = cellfun( @(inds) numel( inds ), burstSpikeInds );
  numSpikesBetweenBursts = cellfun( @(inds) numel( inds ), nonBurstSpikeInds );
  inBurstSpikeRates = numSpikesPerBurst ./ burstDurations;
  betweenBurstSpikeRates = numSpikesBetweenBursts ./ allNonBurstIntervals;
  
  
  bursts = struct( ...
    'burstSpikeInds', { burstSpikeInds }, ...
    'nonBurstSpikeInds', { nonBurstSpikeInds }, ...
    'burstDurations', burstDurations, ...
    'interBurstIntervals', interBurstIntervals, ...
    'burstPeriods', burstPeriods, ...
    'burstRates', burstRates, ...
    'burstDutyCycle', burstDutyCycle, ...
    'numSpikesPerBurst', numSpikesPerBurst, ...
    'numSpikesBetweenBursts', numSpikesBetweenBursts, ...
    'inBurstSpikeRates', inBurstSpikeRates, ...
    'betweenBurstSpikeRates', betweenBurstSpikeRates ...
  );
end

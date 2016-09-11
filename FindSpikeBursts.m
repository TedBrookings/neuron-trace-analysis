%rate = FindSpikeBursts(dT, v, spikes, varargin)
% Currently only makes a smoothed estimate of spike rate
function burstInfo = FindSpikeBursts(dT, v, spikes, varargin)
  parser = inputParser();
  parser.addParameter('filterScale', 3);
  parser.addParameter('plot', true)
  parser.addParameter('minNumSpikes', 10)
  parser.addParameter('maxFilterWidth', 500) % ms
  parser.addParameter('transition', '', ...
                      @(x) ismember( x, {'isi', 'rate', 'times', '', []} ) )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  if ~exist( 'spikes', 'var' ) || isempty( spikes )
    spikes = GetSpikes( dT, v, 'plotSubject', options.plot, ...
                        'debugPlots', options.plot );
  end

  if numel( spikes.times ) < options.minNumSpikes
    spikeStartInd = []; spikeStopInd = [];
  else
    twoPopIsi = fitIsiCutoffTwoPopulation( spikes.times, options )
    twoPopRate = fitRateCutoffTwoPopulation( spikes.times, options )
    twoPopTimes = fitTimesCutoffTwoPopulation( dT, v, spikes, options )
    twoPopZeroDensity = zeroDensityCutoffTwoPopulation( dT, v, spikes, options )
    
    switch options.transition
      case 'isi'
        twoPop = twoPopIsi;
      case 'rate'
        twoPop = twoPopRate;
      case 'times'
        twoPop = twoPopTimes;
      otherwise
        correct = [twoPopIsi.fractionCorrect, ...
                   twoPopRate.fractionCorrect, ...
                   twoPopTimes.fractionCorrect ];
        [maxCorrect, maxInd] = max( correct, [], 'omitnan' );
        switch maxInd
          case 1, twoPop = twoPopIsi;
          case 2, twoPop = twoPopRate;
          case 3, twoPop = twoPopTimes;
        end
        if maxCorrect < 0.8
          twoPop.isiCutoff = -Inf;
          twoPop.isiCutoffFast = -Inf;
          twoPop.isiCutoffSlow = Inf;
        end
    end
    if isfinite( twoPopZeroDensity.isiCutoff ) ...
        && twoPopZeroDensity.isiCutoff < twoPop.isiCutoff
      twoPop = twoPopZeroDensity;
    end
  
    isi = diff( spikes.times );
    [spikeStartInd, spikeStopInd] = ...
      getBurstInds( isi, twoPop, options );
  end
  
  startInd = spikes.n1List(spikeStartInd);
  stopInd = spikes.n2List(spikeStopInd);
  startTime = dT .* (startInd - 1);
  stopTime = dT .* (stopInd - 1);
  burstInfo = struct( 'spikeStartInd', spikeStartInd, ...
                      'spikeStopInd', spikeStopInd, ...
                      'startInd', startInd, ...
                      'stopInd', stopInd, ...
                      'startTime', startTime, ...
                      'stopTime', stopTime );
  
  if options.plot
    rate = getSpikeRate( dT, v, spikes, options );
  
    t = (dT/1000) .* (0:numel(v)-1);
    fig = NamedFigure('SpikeRate'); fig.WindowStyle = 'docked'; clf(fig)
    ax1 = subplot(2,5,1:4, 'Parent', fig);
    plot(ax1, t, v)
    hold( ax1, 'on' )
    minV = min( v ); maxV = max( v ); h = maxV - minV;
    for n = 1:numel( spikeStartInd )
      t1 = startTime(n) / 1000; t2 = stopTime(n) / 1000;
      rectangle( 'Position', [t1, minV, t2-t1, h], 'EdgeColor', 'r', ...
                 'FaceColor', [1 0 0 0.2], 'Parent', ax1 )
    end
    axis( ax1, 'tight' )
    ylabel('voltage (mV)')

    ax2 = subplot(2,5,6:9, 'Parent', fig);
    plot(ax2, t, rate)
    axis( ax2, 'tight' )
    linkaxes( [ax1, ax2], 'x' )
    ylabel('spike rate (Hz)')
    xlabel('time (sec)')
    
    ax3 = subplot(2,5,5:5:10, 'Parent', fig);
    numBins = min( 100, round( sqrt( numel( rate ) ) ) );
    h = histogram( ax3, rate, numBins, 'Normalization', 'probability' );
    axis( ax3, 'tight' )
    ylim( ax3, [0, max( h.Values(2:end) )] );    
  end
  
  if nargout == 0
    varargout = {};
  else
    varargout = {burstInfo};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = fitIsiCutoffTwoPopulation( times, options )
  isi = diff( times );

  [isiDensity, isis] = KernelDensity( isi, [], 0.4 );
  isiDensity = isiDensity ./ sum( isiDensity );
  sigmaIsis = std( isi );
  isi1 = isis(1) + sigmaIsis;
  [~, maxInd] = max( isiDensity );
  if isis(maxInd) < isi1
    isi1 = isis(maxInd);
  end
  isi2 = isis(end) - sigmaIsis;
  while isi2 <= isi1
    isi1 = (isi1 + isis(1)) / 2; isi2 = (isi2 + isis(end)) / 2;
  end
  sigma1 = sigmaIsis; sigma2 = sigmaIsis;
  p1 = sum( isiDensity(1:round(numel(isis)/2)) );
  startParams = [ p1, isi1, sigma1, isi2, sigma2 ];
  
  fitErr = @(params) modelError( params, isis, isiDensity );
  fitParams = fminsearch( fitErr, startParams );
  fprintf( 'fitParams =' ); fprintf( ' %f', fitParams );
  fprintf( '\n' )

  % find dividing line where bursting / non-bursting are equally likely
  isiCutoff = solveCutoff( fitParams, sigmaIsis, 1.0 );
  % find dividing line where bursting is 10 times more likely
  isiCutoffFast = solveCutoff( fitParams, sigmaIsis, 10.0 );
  % find dividing line where non-bursting is 10 times more likely
  isiCutoffSlow = solveCutoff( fitParams, sigmaIsis, 0.10 );
  
  p1 = fitParams(1); r1 = fitParams(2); s1 = fitParams(3);
  p2 = 1.0 - p1; r2 = fitParams(4); s2 = fitParams(5);
  mix1 = p1 .* exp( -((isis-r1)/s1).^2 );
  mix2 = p2 .* exp( -((isis-r2)/s2).^2 );
  mixTot = sum( mix1 + mix2 );
  mix1 = mix1 ./ mixTot; mix2 = mix2 ./ mixTot;

  correct = estimateCorrect( fitParams, isiCutoff );

  twoPop = struct( ...
    'isiCutoff', isiCutoff, ...
    'isiCutoffFast', isiCutoffFast, ...
    'isiCutoffSlow', isiCutoffSlow, ...
    'isis', isis, ...
    'isiDensity', isiDensity, ...
    'fitParams', fitParams, ...
    'mix1', mix1, ...
    'mix2', mix2, ...
    'fractionCorrect', correct ...
  );

  if options.plot
    fig = NamedFigure( 'Isi Kernel Burst' ); fig.WindowStyle = 'docked';
    clf( fig )
    ax = subplot( 1,1,1, 'Parent', fig ); hold( ax, 'on' )
    plot( ax, twoPop.isis, twoPop.isiDensity, 'b-' );
    plot( ax, twoPop.isis, finiteMixtureModel( twoPop.fitParams, twoPop.isis ), 'r-' )
    plot( ax, [isiCutoff isiCutoff], [0 max(twoPop.isiDensity)], 'k--' )
    plot( ax, twoPop.isis, twoPop.mix1, 'm:' )
    plot( ax, twoPop.isis, twoPop.mix2, 'm:' )
    legend( ax, {'ISIs', 'Model fit', 'Cutoff', 'Populations'}, ...
            'Location', 'Best' )
  end

  varargout = {twoPop};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function correct = estimateCorrect( fitParams, cutoff )
  if ~isreal( cutoff ) || ~isfinite( cutoff )
    correct = 0;
    return
  end
  p1 = fitParams(1); center1 = fitParams(2); s1 = fitParams(3);
  p2 = 1.0 - p1; center2 = fitParams(4); s2 = fitParams(5);
  correct1 = (p1 / 2) * (1 + erf( (cutoff - center1 ) / s1 / sqrt( 2 ) ));
  correct2 = (p2 / 2) * (1 - erf( (cutoff - center2 ) / s2 / sqrt( 2 ) ));
  correct = correct1 + correct2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = fitRateCutoffTwoPopulation( times, options )
  isi = diff( times );
  [isiDensity, isis] = KernelDensity( isi, [], 0.4 );
  rate = 1.0 ./ isi; rates = linspace( 0, max( rate ), 1000 );
  [rateDensity, rates] = KernelDensity( rate, rates, 0.4 );
  rateDensity = rateDensity ./ sum( rateDensity );
  
  sigmaRates = std( rates );
  rate2 = rates(end) - sigmaRates;
  [~, maxInd] = max( rateDensity );
  if rates(maxInd) > rate2
    rate2 = rates(maxInd);
  end
  rate1 = rates(1) + sigmaRates;
  while rate2 <= rate1
    rate1 = (rate1 + rates(1)) / 2; rate2 = (rate2 + rates(end)) / 2;
  end
  halfInd = round( numel( rates ) / 2 );
  sigma1 = sigmaRates; sigma2 = sigmaRates;
  p1 = sum( rateDensity(1:halfInd) );
  startParams = [ p1, rate1, sigma1, rate2, sigma2 ];
  
  fitErr = @(params) modelError( params, rates, rateDensity );
  fitParams = fminsearch( fitErr, startParams );
  fprintf( 'fitParams =' ); fprintf( ' %f', fitParams );
  fprintf( '\n' )

  % find dividing line where bursting / non-bursting are equally likely
  rateCutoff = solveCutoff( fitParams, sigmaRates, 1.0 );
  isiCutoff = 1.0 / rateCutoff;

  % find dividing line where bursting is 10 times more likely
  isiCutoffFast = 1.0 / solveCutoff( fitParams, sigmaRates, 0.05 );
  % find dividing line where non-bursting is 10 times more likely
  isiCutoffSlow = 1.0 / solveCutoff( fitParams, sigmaRates, 20.0 );
  
  p1 = fitParams(1); r1 = fitParams(2); s1 = fitParams(3);
  p2 = 1.0 - p1; r2 = fitParams(4); s2 = fitParams(5);
  mix1 = p1 .* exp( -((rates-r1)/s1).^2 );
  mix2 = p2 .* exp( -((rates-r2)/s2).^2 );
  mixTot = sum( mix1 + mix2 );
  mix1 = mix1 ./ mixTot; mix2 = mix2 ./ mixTot;

  rateCutoff
  p1, r1, s1
  p2, r2, s2
  
  correct = estimateCorrect( fitParams, rateCutoff );
  
  twoPop = struct( ...
    'isiCutoff', isiCutoff, ...
    'isiCutoffFast', isiCutoffFast, ...
    'isiCutoffSlow', isiCutoffSlow, ...
    'isis', isis, ...
    'isiDensity', isiDensity, ...
    'fitParams', fitParams, ...
    'mix1', mix1, ...
    'mix2', mix2, ...
    'fractionCorrect', correct ...
  );

  if options.plot
    fig = NamedFigure( 'Rate Kernel Burst' ); fig.WindowStyle = 'docked';
    clf( fig )
    ax = subplot( 1,1,1, 'Parent', fig ); hold( ax, 'on' )
    plot( ax, rates, rateDensity, 'b-' );
    plot( ax, rates, finiteMixtureModel( twoPop.fitParams, rates ), 'r-' )
    plot( ax, [rateCutoff rateCutoff], [0 max(rateDensity)], 'k--' )
    plot( ax, rates, twoPop.mix1, 'm:' )
    plot( ax, rates, twoPop.mix2, 'm:' )
    legend( ax, {'ISIs', 'Model fit', 'Cutoff', 'Populations'}, ...
            'Location', 'Best' )
  end

  varargout = {twoPop};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = fitTimesCutoffTwoPopulation( dT, v, spikes, options )
  isi = diff( spikes.times );
  [isiDensity, isis] = KernelDensity( isi, [], 0.4 );
  rate = getSpikeRate( dT, v, spikes, options );
  %[rate, t] = KernelDensity( times, [], 0.2 );
  [rateDensity, rates] = KernelDensity( rate, [], 0.1 );
  rateDensity = rateDensity ./ sum( rateDensity );
  
  cumulative = cumsum( rateDensity );
  halfInd = find( cumulative >= 0.5, 1 ); halfRate = rates(halfInd);
  [~, maxInd] = max( rateDensity(1:halfInd) );
  rate1 = rates(maxInd);
  sigma1 = std( rate(rate <= halfRate) );
  
  [~, maxInd] = max( rateDensity(halfInd+1:end) );
  rate2 = rates(maxInd+halfInd);
  sigma2 = std( rate(rate >= halfRate) );
  sigmaRates = std( rates );
  %rate2 = rates(end) - sigmaRates;
  %[~, maxInd] = max( rateDensity );
  %if rates(maxInd) > rate2
  %  rate2 = rates(maxInd);
  %end
  %rate1 = rates(1) + sigmaRates;
  while rate2 <= rate1
    rate1 = (rate1 + rates(1)) / 2; rate2 = (rate2 + rates(end)) / 2;
  end
  halfInd = round( numel( rates ) / 2 );
  %sigma1 = sigmaRates; sigma2 = sigmaRates;
  p1 = sum( rateDensity(1:halfInd) );
  startParams = [ p1, rate1, sigma1, rate2, sigma2 ];
  
  fitErr = @(params) modelError( params, rates, rateDensity );
  fitParams = fminsearch( fitErr, startParams );
  fprintf( 'fitParams =' ); fprintf( ' %f', fitParams );
  fprintf( '\n' )

  % find dividing line where bursting / non-bursting are equally likely
  rateCutoff = solveCutoff( fitParams, sigmaRates, 1.0 );
  isiCutoff = 1000.0 / rateCutoff;

  
  p1 = fitParams(1); r1 = fitParams(2); s1 = fitParams(3);
  p2 = 1.0 - p1; r2 = fitParams(4); s2 = fitParams(5);
  mix1 = p1 .* exp( -((rates-r1)/s1).^2 );
  mix2 = p2 .* exp( -((rates-r2)/s2).^2 );
  mixTot = sum( mix1 + mix2 );
  mix1 = mix1 ./ mixTot; mix2 = mix2 ./ mixTot;

  % find dividing line where bursting is 10 times more likely
  isiCutoffFast = 1000.0 / solveCutoff( fitParams, sigmaRates, 0.05 );
  % find dividing line where non-bursting is 10 times more likely
  isiCutoffSlow = 1000.0 / solveCutoff( fitParams, sigmaRates, 20.0 );
    
  rateCutoff
  p1, r1, s1
  p2, r2, s2
  
  correct = estimateCorrect( fitParams, rateCutoff );
  
  twoPop = struct( ...
    'isiCutoff', isiCutoff, ...
    'isiCutoffFast', isiCutoffFast, ...
    'isiCutoffSlow', isiCutoffSlow, ...
    'isis', isis, ...
    'isiDensity', isiDensity, ...
    'fitParams', fitParams, ...
    'mix1', mix1, ...
    'mix2', mix2, ...
    'fractionCorrect', correct ...
  );

  if options.plot
    fig = NamedFigure( 'Rate vs times Kernel Burst' );
    fig.WindowStyle = 'docked'; clf( fig )
    ax = subplot( 1,1,1, 'Parent', fig ); hold( ax, 'on' )
    plot( ax, rates, rateDensity, 'b-' );
    plot( ax, rates, finiteMixtureModel( twoPop.fitParams, rates ), 'r-' )
    plot( ax, [rateCutoff rateCutoff], [0 max(rateDensity)], 'k--' )
    plot( ax, rates, twoPop.mix1, 'm:' )
    plot( ax, rates, twoPop.mix2, 'm:' )
    legend( ax, {'ISIs', 'Model fit', 'Cutoff', 'Populations'}, ...
            'Location', 'Best' )
  end


  varargout = {twoPop};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = zeroDensityCutoffTwoPopulation( dT, v, spikes, options )
  isi = diff( spikes.times );
  [isiDensity, isis] = KernelDensity( isi, [], 0.4 );

  mix1 = isiDensity; mix2 = isiDensity;
  
  zeroLevel = 0.5 / numel( isi );
  posInd = find( isiDensity > zeroLevel, 1 );
  if isempty( posInd ), posInd = 1; end
  zeroInd = posInd + find( isiDensity(posInd+1:end) < zeroLevel, 1 );
  if ~isempty( zeroInd )
    isiCutoff = isis(zeroInd);
    isiCutoffFast = isiCutoff;
    isiCutoffSlow = isiCutoff;
    mix1(isis < isiCutoff) = 0;
    mix2(isis > isiCutoff) = 0;
  else
    isiCutoff = -Inf;
    isiCutoffFast = -Inf;
    isiCutoffSlow = Inf;
  end
  
   twoPop = struct( ...
    'isiCutoff', isiCutoff, ...
    'isiCutoffFast', isiCutoffFast, ...
    'isiCutoffSlow', isiCutoffSlow, ...
    'isis', isis, ...
    'isiDensity', isiDensity, ...
    'fitParams', [], ...
    'mix1', mix1, ...
    'mix2', mix2, ...
    'fractionCorrect', NaN ...
  );
  varargout = {twoPop};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = modelError( params, rates, rateDensity )
  err = sum( (rateDensity - finiteMixtureModel( params, rates )).^2 );
  if params(1) < 0
    err = err + params(1)^2;
  elseif params(1) > 1
    err = err + (params(1) - 1)^2;
  end
  if params(2) < -abs( params(3) )
    err = err + (params(2) + abs( params(3) ))^2;
  end
  if params(3) < 0
    err = err + params(3)^2;
  end
  if params(4) < params(2)
    err = err + (params(4) - params(2))^2;
  end
  if params(5) < 0
    err = err + params(5)^2;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rateDensity = finiteMixtureModel( params, rates )
  % params = [ prob1, rate1, sigma1, rate2, sigma2 ]
  prob1 = params(1); rate1 = params(2); sigma1 = params(3);
  if rate1 < 0
    rate1 = 0;
  end
  if prob1 < 0
    prob1 = 0;
  elseif prob1 > 1
    prob1 = 1;
  end
  prob2 = 1.0 - prob1; rate2 = params(4); sigma2 = params(5);
  rateDensity = prob1 .* exp( -((rates - rate1)./sigma1).^2 ) ...
              + prob2 .* exp( -((rates - rate2)./sigma2).^2 );
  rateDensity = rateDensity ./ sum( rateDensity );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p1 exp( -(r-r1)^2/s1^2 ) ~ ratio * p2 exp( -(r-r2)^2/s2^2 )
function isiCutoff = solveCutoff( fitParams, sigmaIsis, ratio )
  % p1 exp( -(r-r1)^2/s1^2 ) ~ ratio * p2 exp( -(r-r2)^2/s2^2 )
  % -log( p2 * ratio / p1 ) = (r - r1)^2/s1^2 - (r-r2)^2/s2^2
  % r^2 [1/s1^2 - 1/s2^2] + r*2*[r2/s2^2 - r1/s1^2] + r1^2/s1^2 - r2^2/s2^2 + log(p2 * ratio/p1) = 0
  p1 = fitParams(1); r1 = fitParams(2); s1 = fitParams(3);
  p2 = 1.0 - p1; r2 = fitParams(4); s2 = fitParams(5);
  a = 1 / s1^2 - 1 / s2^2;
  b = 2 * (r2 / s2^2 - r1 / s1^2);
  c = r1^2 / s1^2 - r2^2 / s2^2 + log( ratio * p2 / p1 );
  if abs( a * (s1^2 + s2^2) ) > 1e-3
    root = b^2 - 4*a*c;
    if root > 0
      root = sqrt( root );
      isiCutoff1 = (-b + root) / (2 * a );
      isiCutoff2 = (-b - root) / (2 * a );
    
      if isiCutoff1 > 0 && isiCutoff2 < 0
        isiCutoff = isiCutoff1;
      elseif isiCutoff1 < 0 && isiCutoff2 > 0
        isiCutoff = isiCutoff2;
      elseif r1 <= isiCutoff1 && r2 >= isiCutoff1
        isiCutoff = isiCutoff1;
      else
        isiCutoff = isiCutoff2;
      end
    else
      isiCutoff = -b / (2 * a );
    end
  else
    isiCutoff = -c / b;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [spikeStartInd, spikeStopInd] = ...
    getBurstInds( isi, twoPop, options )
  numIsi = numel( isi );

  inBurst = false;
  spikeStartInd = []; spikeStopInd = []; stopInd = 0;
  for n = 1:numIsi
    if inBurst
      if isi(n) > twoPop.isiCutoffSlow
        startInd = spikeStartInd(end);
        stopInd = startInd - 1 ...
                + find( isi(startInd:n) <= twoPop.isiCutoff, 1, 'last' );
        spikeStopInd = [ spikeStopInd, 1 + stopInd ]; %#ok<AGROW>
        inBurst = false;
        stopInd = n;
      end
    else
      if isi(n) > twoPop.isiCutoffSlow
        stopInd = n;
      elseif isi(n) < twoPop.isiCutoffFast
        startInd = stopInd ...
                 + find( isi(stopInd+1:n) <= twoPop.isiCutoff, 1 );
        spikeStartInd = [ spikeStartInd, startInd ]; %#ok<AGROW>
        inBurst = true;
      end
    end
  end
  
  if inBurst
    spikeStopInd = [spikeStopInd, numIsi];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rate = getSpikeRate( dT, v, spikes, options )
  rate = zeros(size(v));
  if numel( spikes.n1List ) < options.minNumSpikes
    return
  end
  rate(spikes.maxV.ind) = 1000.0/dT;
  
  filtW = options.filterScale * median( diff( spikes.n1List ) );
  filtW = min( filtW, options.maxFilterWidth / dT );
  halfFiltLen = 3 * round( filtW );
  
  filt = exp( -(-halfFiltLen:halfFiltLen).^2 ./ (2 * filtW^2) );
  filt = filt ./ sum( filt );
  
  rate = applyFilter( rate, filt );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply convolution filter
function y = applyFilter( y, filt )
  transpose = ~isrow( y );
  if transpose
    y = y';
  end
  filtLen = numel(filt);
  halfLen = (filtLen - 1) / 2;
  % pad y symmetrically
  y = [ y(halfLen+1:-1:2), y, y(end-1:-1:end-halfLen) ];
  % apply convolution filter, keeping valid part of y
  y = conv( y, filt, 'valid' );
  if transpose
    y = y';
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diffy = applyFiltNoDelay(y, dFilt)
% apply convolution filter

% first make y periodic
fLen = length(dFilt.filter);
nHalf = (fLen - 1) / 2;
y1 = y(1); yN = y(end); yNonperiodic = linspace(y1, yN, length(y));
y = y - yNonperiodic;
% next pad y with appropriate symmetry
if mod(dFilt.order, 2) == 0
  % order is even, so pad symmetrically
  y = [y(nHalf+1:-1:2), y, y(end-1:-1:end-nHalf)];
else
  % order is odd, so pad antisymmetrically
  y = [-y(nHalf+1:-1:2), y, -y(end-1:-1:end-nHalf)];
end

% apply convolution filter
diffy = conv(y, dFilt.filter, 'valid');
% apply any needed corrections resulting from making y periodic
if dFilt.order == 1
  averageSlope = (yN - y1) / ( dFilt.dx * (length(y) - 1) );
  diffy = diffy + averageSlope;
elseif dFilt.order == 0
  diffy = diffy + yNonperiodic;
end
end

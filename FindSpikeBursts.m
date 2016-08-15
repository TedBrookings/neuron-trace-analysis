%rate = FindSpikeBursts(dT, v, spikes, varargin)
% Currently only makes a smoothed estimate of spike rate
function varargout = FindSpikeBursts(dT, v, spikes, varargin)
  parser = inputParser();
  parser.addParameter('filterScale', 3);
  parser.addParameter('plot', true)
  parser.addParameter('minNumSpikes', 10)
  parser.addParameter('maxFilterWidth', 500) % ms
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  if ~exist('spikes', 'var') || isempty(spikes)
    spikes = GetSpikes( dT, v, 'plotSubject', options.plot, ...
                        'debugPlots', options.plot );
  end

  
  rate = getSpikeRate( dT, v, spikes, options );

  if options.plot
    t = (dT/1000) .* (0:numel(v)-1);
    fig = NamedFigure('SpikeRate'); fig.WindowStyle = 'docked'; clf(fig)
    ax1 = subplot(2,5,1:4, 'Parent', fig);
    plot(ax1, t, v)
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
    varargout = {rate};
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
  if ~isrow(y)
    transpose = true;
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

function varargout = GetAverageSpikeWaveform( trace, n1List, n2List, ...
                                              varargin )
  parser = inputParser();
  parser.addParameter( 'plot', false )
  parser.addParameter( 'title', 'Average Spike Waveform' )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  if ~isrow( trace ), trace = trace'; end
  allWaveforms = alignWaveforms( trace, n1List, n2List );

  meanWaveform = mean( allWaveforms, 1, 'omitnan' );
  varargout = { meanWaveform, allWaveforms };
  
  if options.plot
    w = numel( meanWaveform ); wHalf = (w-1)/2;
    x = -wHalf:wHalf;
    fig = NamedFigure( options.title ); clf( fig );
    fig.WindowStyle = 'docked';
    ax = axes( 'Parent', fig, 'LooseInset', [0 0 0 0], ...
               'OuterPosition', [0 0 1 1] );
    hold( ax, 'on' )
    plot( ax, x, allWaveforms' );
    plot( ax, x, meanWaveform, 'k:', 'LineWidth', 4 )
    ylabel( ax, 'Voltage' );
    xlabel( ax, 'sample number' )
    axis( ax, 'tight' )
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align spike waveforms into numSpikes x w matrix with spike peak at center 
function allWaveforms = alignWaveforms( trace, n1List, n2List )
  w = median( n2List - n1List );
  wHalf = ceil( (w-1) / 2 ); mid = wHalf + 1; w = 2 * wHalf + 1;
  numSpikes = numel( n1List );
  for ind = numSpikes:-1:1
    n1 = n1List(ind); n2 = n2List(ind);
    allWaveforms(ind,:) = alignSpike( trace(n1:n2), mid, w );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align this spike so that spike peak is at center
function waveform = alignSpike( waveform, mid, w )
  [~, maxInd] = max( waveform );
  delta = maxInd - mid;
  if delta > 0
    waveform(1:delta) = [];
  elseif delta < 0
    waveform = [ NaN( 1, -delta), waveform ];
  end
  delta = numel( waveform ) - w;
  if delta > 0
    waveform(w+1:end) = [];
  else
    waveform = [waveform, NaN( 1, -delta )];
  end
end

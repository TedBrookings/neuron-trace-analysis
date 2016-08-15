function trace = RemoveSpikes(dt, trace, varargin)
% trace = RemoveSpikes(trace, varargin)
% Remove spikes from trace by interpolating between beginning and end of
% spikes

parser = inputParser();
parser.KeepUnmatched = true;
parser.StructExpand = false;
parser.addParameter('spikes', [])
parser.addParameter('spikeTrace', [])
parser.addParameter('plot', false)

parser.parse( varargin{:} )
options = parser.Results;
unmatchedFields = fieldnames( parser.Unmatched );
for n = 1:numel(unmatchedFields)
  options.(unmatchedFields{n}) = parser.Unmatched.(unmatchedFields{n});  
end

if isempty(options.spikes)
  if isempty(options.spikeTrace)
    spikeTrace = trace;
  else
    spikeTrace = options.spikeTrace;
  end
  spikes = GetSpikes( dt, spikeTrace, options );
else
  spikes = options.spikes;
end


if options.plot
  fig = NamedFigure('RemoveSpikes') ; fig.WindowStyle = 'docked'; clf(fig)
  t = (0.001 * dt) * (0:(length(trace)-1));
  ax = subplot(1,1,1, 'Parent', fig );
  plot( ax, t, trace )
  trace = removeSpikesFromTrace( trace, dt, spikes );
  hold( ax, 'on' )
  plot( ax, t, trace, 'r-' )
  hold( ax, 'off' )
  axis( ax, 'tight' )
else
  trace = removeSpikesFromTrace( trace, dt, spikes );
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trace = removeSpikesFromTrace( trace, dt, spikes )
%n1Barrier = 1;
maxSpacing = 200.0 / dt;
for m = 1:length(spikes.n1List)
  n1A = spikes.n1List(m); n2A = spikes.n2List(m);
  iA = spikes.preMaxCurve.ind(m);
  if iA > 1
    iB = iA ; iA = iA - 1;
  else
    iB = iA + 1;
  end

  if m < length(spikes.n1List) && ...
      (spikes.n1List(m+1) - n2A) < maxSpacing
    % join from beginning of this spike to beginning of next spike
    iD = spikes.preMaxCurve.ind(m+1);
    iC = iD - 1;
    interpType = 'linear';
  else
    % join from beginning of this spike to end of this spike
    iD = spikes.postMinV.ind(m);
    iC = iD - 1;
    interpType = 'linear';
  end
  x = [iA iB iC iD]; y = trace(x); xInterp = (iB+1):(iC-1);
  trace(xInterp) = interp1(x, y, xInterp, interpType);
  %n1Barrier = n2B;
end
end
% fig = NamedFigure( name, varargin )
%  Looks for a figure with requested name, sets it to the current figure,
%  and returns a handle fig to the figure.  If no such figure exists, it
%  creates one, names it, and returns a handle
%    INPUT:
%     -name   The name of the figure
%    OPTIONS:
%      additional key,value pairs will be passed to the figure() function,
%      allowing the setting of figure properties at creation time
%    OUTPUT:
%     -fig    Handle to the figure
function varargout = NamedFigure( name, varargin )
  parser = inputParser();
  parser.KeepUnmatched = true;
  parser.addParameter( 'Position', [] )
  parser.addParameter( 'WindowStyle', 'normal' )
  parser.addParameter( 'Visible', 'on' )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  % find any figures with requested name
  fig = findobj( 'Name', name, 'Type', 'Figure' );
  
  defaultProperties = { ...
    'PaperSize', [10 7.5], ...
    'InvertHardcopy', 'off' ...
  };
  if isempty( fig )
    % no such figure exists, create one
    
    % set some default properties for figures
    set( 0, 'defaultaxesfontname', 'Arial' )
    set( 0, 'defaulttextfontname', 'Arial' )
    set( 0, 'defaultaxesfontsize', 15 )
    set( 0, 'defaulttextfontsize', 18 )
    
    % create the figure
    fig = figure( 'Name', name, defaultProperties{:}, varargin{:}, ...
                  'visible', 'off', 'WindowStyle', 'normal' );
    if isempty( options.Position )
      units = fig.Units;
      fig.Units = 'inches';
      fig.Position = [0, 0, fig.PaperSize];
      fig.Units = units;
    end
    if strcmp( options.WindowStyle, 'docked' )
      fig.WindowStyle = 'docked';
    end
    if strcmp( options.Visible, 'on' )
      fig.Visible = 'on';
    end
  else
    % figure exists, set it to current figure
    set( 0, 'CurrentFigure', fig(1) )
    if nargin > 1
      % update figure properties
      assert( mod( numel( varargin ), 2 ) == 0, ...
              'Figure options should be key, value pairs' )
      for n = 1:2:numel( varargin )
        fig.(varargin{n}) = varargin{n+1};
      end
    end
  end

  % return the figure handle if it's requested
  if nargout == 0
    varargout = {};
  else
    varargout = { fig };
  end
end
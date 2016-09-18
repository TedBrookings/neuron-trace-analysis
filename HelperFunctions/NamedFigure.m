% fig = NamedFigure( name, varargin )
%  Looks for a figure with requested name, sets it to the current figure,
%  and returns a handle fig to the figure.  If no such figure exists, it
%  creates one, names it, and returns a handle
%    INPUT:
%     -name   The name of the figure
%    OUTPUT:
%     -fig    Handle to the figure
function varargout = NamedFigure( name, varargin )
  % find any figures with requested name
  fig = findobj( 'Name', name, 'Type', 'Figure' );
  
  if isempty( fig )
    % no such figure exists, create one
    
    % set some default properties for figures
    set( 0, 'defaultaxesfontname', 'Arial' )
    set( 0, 'defaulttextfontname', 'Arial' )
    set( 0, 'defaultaxesfontsize', 15 )
    set( 0, 'defaulttextfontsize', 18 )
    
    % create the figure
    fig = figure( 'Name', name, varargin{:} );
    
    % set convenient figure properties
    set( fig, 'Name', name )
  else
    % figure exists, get it
    set( 0, 'CurrentFigure', fig )
    if nargin > 1
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
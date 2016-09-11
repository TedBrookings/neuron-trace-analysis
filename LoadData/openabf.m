% This file tells MATLAB how to open .abf files when they are dragged onto
% the main workspace
function output = openabf( fileName )
  output = LoadAbf( fileName );
  assignin( 'base', 'abfData', output )
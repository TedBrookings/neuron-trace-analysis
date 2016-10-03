% homeStr = Home()
% return string with path to home directory of user
function homeStr = Home()
  if ispc
    homeStr = [ getenv( 'HOMEDRIVE' ), getenv( 'HOMEPATH' ) ];
  else
    homeStr = getenv( 'HOME' );
  end
end
% Add binica to Matlab's path
path_startup     = userpath(); 
path_startup     = path_startup(1:end-1);
path_binica      = fullfile(path_startup, 'binica_osx_fat');
path_environment = getenv('PATH');
path_environment = [path_environment ':' path_binica];
setenv('PATH', path_environment);
[status, result] = system('echo $PATH');


% Test if Binica is installed
[status, result] = system('ica_osx');
if(isempty(strfind(result, 'command not found')))
    display('        Binica Installed');
else
    display('Binica Not Installed');
end
% Add binica to Matlab's path
[path_eeglab, ~] = fileparts(which('eeglab'));
path_binica      = fullfile(path_eeglab, 'binica_osx_fat');
path_environment = getenv('PATH');
path_environment = [path_environment ':' path_binica];

% Add BINICA to Unix path
setenv('PATH', path_environment);
% [status, result] = system('echo $PATH');

% Add BINICA to Matlab path
addpath(path_binica);


% Test if Binica is installed
[status, result] = system('ica_osx');
if(strfind(result, 'command not found'))
    display('        Binica not installed');
elseif(strfind(result, 'Permission denied'))
    display('        Binica not installed: Permission denied');
else
    
    if(exist('ica_osx') == 2) %#ok<EXIST>
        display('Binica installed');
    else
        display('Binica not installed: INSERT REASON');
    end
end
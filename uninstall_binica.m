

%% Remove `binica_osx_fat` folder
foldername_binica   = 'binica_osx_fat';
[path_eeglab, ~]    = fileparts(which('eeglab'));
path_binica         = fullfile(path_eeglab, foldername_binica);
system(['rm -rf ', path_binica])

%% Remove startup file
path_matlab      = userpath();           % Get Matlab's startup directory
path_matlab      = path_matlab(1:end-1); % remove colon (:) at end of string
filename_startup = 'startup.m';

path_startup = fullfile(path_matlab, filename_startup);
if (exist(path_startup, 'file'))
    % Append startup info
    system(['mv ', path_startup, ' ', path_startup, '.old'])
else
    % Copy startup.m file to Matlab startup directory
    system(['cp ', filename_startup, ' ', path_matlab]);
end




% Unzip 'Binica_mac_intel.zip' to eeglab directory
binica_zip = 'Binica_mac_intel.zip';
[path_eeglab, ~] = fileparts(which('eeglab'));
system(['unzip ' binica_zip , ' -d ' , path_eeglab])


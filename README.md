# Install BINICA

To install:

- **In Matlab, run `install_binica`**
  - This will unzip the `binica_osx_fat` folder into your EEGLAB folder, and add it to your Matlab path and Unix path
- **Edit `icadefs.m`**:  Set the `ICABINARY` variable to `'ica_osx'` binary file.
- `.../eeglab13_x_x/functions/sigprocfunc/icadefs.m` (line 119)
  - ```matlab
    ...
    % INSERT location of ica executable (LINUX ONLY) for binica.m below
    eeglab_p = fileparts(which('eeglab'));
    ICABINARY = 'ica_osx';  % <<< EDIT HERE
    
    try
       set(0,'defaultaxesfontsize',AXES_FONTSIZE);
    ...

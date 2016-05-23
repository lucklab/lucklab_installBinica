# Install BINICA

To install:

Run LuckLab BINICA Installer:
- Open Matlab
- Change to `lucklab_installBinica` directory
- Run the `install_binica` script


Edit `icadefs.m` file
- Edit the `icadefs.m file` to point to the `ica_osx` executable you intend to use (variable ICABINARY).
- `.../eeglab13_x_x/functions/sigprocfunc/icadefs.m` (line 119)
  - ```matlab
    ...
    % INSERT location of ica executable (LINUX ONLY) for binica.m below
    eeglab_p = fileparts(which('eeglab'));
    ICABINARY = 'ica_osx'; 
    
    try
       set(0,'defaultaxesfontsize',AXES_FONTSIZE);
    ...

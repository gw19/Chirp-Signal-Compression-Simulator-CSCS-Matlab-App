# CSCS
### CSCS (Chirp Signal Compression Simulator) is a Matlab application for easily simulating the single/multi-target linear frequency modulation (chirp) signal compression.

https://www.mathworks.com/matlabcentral/fileexchange/58719-chirp-signal-compression-simulator-cscs

#### Double click the .mlappinstall file to install the Matlab application.
#### Need Matlab R2016a or newer version.
#### Without any Matlab toolbox.

Key words: signal, matched filter, compression, linear frequency modulation (linear FM), chirp, RADAR.  

**Step 1:** Chirp Set - Set the basic chirp.  
First, determine the bandwidth, pulse duration, sampling frequency and chirp form you want, it will display a simple basic chirp waveform s(t), and results of the output signal compressed after a default, corresponding matched filter h(t). You can also add White Gaussian Noises into the original signals, and Kaiser Window into the matched filter. It will auto update the results when you change the initial setting. Then it will use this basic waveform for both of the single and multiple targets simulation.  

**Step 2:** Targets - Set ranges of targets.  
Second, assume that all of the targets are on the ground range, input the radar height, incidence angle, and ranges between the near range and targets, check if it can be identified by compressing the combination signal after the matched filter.  

Thank You,  
GW  
imgw19@gmail.com  

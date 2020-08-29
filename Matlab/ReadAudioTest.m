audioinfo = audiodevinfo;

recObj = audiorecorder;
Fs = 44100 ; 
nBits = 16 ; 
nChannels = 1 ; 
ID = 1; % default audio input device 
recObj = audiorecorder(Fs,nBits,nChannels,ID);

disp('Start speaking.')
recordblocking(recObj,5);
disp('End of Recording.');
play(recObj);
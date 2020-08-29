<<<<<<< HEAD
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
=======
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
>>>>>>> 01bbdc2f27f3e182adf89ab74ea82be8eaac7960
play(recObj);
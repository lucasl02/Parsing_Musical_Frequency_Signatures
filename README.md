# Parsing_Musical_Frequency_Signatures
In this project, spectral analysis on an audio file was done to isolate
specific instrument in the track. Specifically, the Gabor Transform was
applied to create a spectrogram of the audio file. A Gaussian filter
centered about the peak frequencies was then used to isolate these
frequencies, which were then further filtered by a make-shift Shannon
filter to single out the bass and guitar portions of the audio.

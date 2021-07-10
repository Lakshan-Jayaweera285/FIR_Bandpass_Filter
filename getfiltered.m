function[out] = getfiltered(xnt,filter)
Npoint = length(xnt) + length(filter) - 1;  % length for fft in x dimension
xfft = fft(xnt,Npoint);
filterfft = fft(filter,Npoint);
outfft = filterfft .* xfft;
out = ifft(outfft,Npoint);

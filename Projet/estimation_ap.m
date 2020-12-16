function tab_ap = estimation_ap (signal,p)
corrSig = xcorr(signal,p);

rxx_vect = fliplr(corrSig(1:p+1));

MatToeplitz = toeplitz(rxx_vect(1:end-1));
tab_ap = -inv(MatToeplitz)*rxx_vect(2:end)';
end
function P=FFT_spectrum_generate(x,N)
   Y=fft(x,N);
   P1=abs(Y/N); 
   P = P1(1:N/2+1);
   P(2:end-1) = 2*P(2:end-1);
end
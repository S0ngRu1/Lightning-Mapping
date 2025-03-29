function data_original=window_plus(N,data_original)
window_hann = hann(N*0.2);
data_original(1:N*0.1)=data_original(1:N*0.1).*window_hann(1:N*0.1);
data_original(N*0.9+1:end)=data_original(N*0.9+1:end).*window_hann(N*0.1+1:end);
end
function data_filter = datafilter(data_original1)
    Wc=2*80e6/200e6;
    Wp=2*30e6/200e6;
    [b,a] = butter(2,[Wp,Wc]);
    data_filter = (filtfilt(b,a,data_original1));
end
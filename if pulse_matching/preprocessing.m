function preprocessed = preprocessing(signal)
    signal1_removed = detrend(signal);
    preprocessed = windowsignal(signal1_removed);
end
signal_length = 1024;
r_loction = 4e8;

ch1 = read_signal('..\\20230718175104.9180CH1.dat',signal_length,r_loction);
ch2 = read_signal('..\\20230718175104.9180CH2.dat',signal_length,r_loction);
ch3 = read_signal('..\\20230718175104.9180CH3.dat',signal_length,r_loction);

upsampling_factor = 50;
[ch1_up, ch2_up, ch3_up] = deal(...
        upsampling(ch1, upsampling_factor,'polyfit'), ...
        upsampling(ch2, upsampling_factor,'polyfit'), ...
        upsampling(ch3, upsampling_factor,'polyfit'));
    ch1_upsp = ch1_up(:,2);
    ch2_upsp = ch2_up(:,2);
    ch3_upsp = ch3_up(:,2);

  plot(ch1)
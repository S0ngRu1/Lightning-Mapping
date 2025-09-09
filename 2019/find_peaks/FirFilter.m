function filtered_signal = FirFilter(signal)
   resp = 'bandpassfir';  %选择滤波器类型
   order = 50;   
   CutoffFrequency1 = 20e6; 
   CutoffFrequency2 = 80e6;
   Fs = 200e6;              %设置信号采样频率
   compFlag = 'on';     %compFlag：是否进行延迟补偿，'on'为开启补偿，'off'为关闭补偿。      
 
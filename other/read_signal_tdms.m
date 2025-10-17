%% %%%tdms 文件读取 引雷点   
 function signal = read_signal_tdms(signal_path, r_length,r_loction)
[output_400,~] = TDMS_getStruct(signal_path);    
int_400=output_400.groups.chans.data; %%%数据
signal = double(int_400(r_loction:r_length+r_loction-1))';
end

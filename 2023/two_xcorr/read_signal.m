function signal = read_signal(signal_path, r_length)
    fid  = fopen(signal_path,'r');
    r_loction = 4401e5;%读取数据的位置

    %使用fseek函数将文件指针移动到指定位置，以便读取数据。
    %这里指定移动位置为r_location，表示移动到指定位置开始读取数据。
    fseek(fid,r_loction*2,'bof');
    %使用fread函数从文件中读取数据，读取的数据长度为r_length，数据以int16格式读取。
    %将读取到的数据分别保存到变量ch_1、ch_2和ch_3中。
    signal = fread(fid,r_length,'int16');
    %关闭所有文件
    fclose('all');
end


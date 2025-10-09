
%% CHJ 干涉仪400 2025
clear
%%%数据存在格式为两个通道交替
%%%如10个数据 1 2 3 4 5 6 7 8 9 10
%%%通道1内数据为 1 3 5 7 9 
%%%通道2内数据为 2 4 6 8 10
%%%数据为int16格式    每个占2个字节
%%%%%%如果需要分段读取，计算好间隔，不然两个通道会混合
skip_byte=2; %%每个数据占2个字节，因此该参数为2时，提取波形为400Mhz。改变为2*n，采样率减小为400/n
%%%CardOne 混合保存通道1和通道2数据
int_filename_CardOne='software_CardOne_20250823172542.dat';
read_int_CardOne=strcat('H:\CHJ400\data\',int_filename_CardOne);
fid_CardOne=fopen(read_int_CardOne,'rb');
where_CardOne=fseek(fid_CardOne,0, 'bof'); 
CH1_data = fread(fid_CardOne, inf, 'int16',skip_byte,'b');   
where_CardOne=fseek(fid_CardOne,2, 'bof'); 
CH2_data = fread(fid_CardOne, inf, 'int16',skip_byte,'b'); 
% CH1_index=1:2:length(data_CardOne_read);
% CH2_index=2:2:length(data_CardOne_read);
% 
% CH1_data=data_CardOne_read(CH1_index);
% CH2_data=data_CardOne_read(CH2_index);
%%%CardTwo 混合保存通道3和通道4数据
int_filename_CardTwo='software_CardTwo_20250823172542.dat';
read_int_CardTwo=strcat('H:\CHJ400\data\',int_filename_CardTwo);
fid_CardTwo=fopen(read_int_CardTwo,'rb');
where_CardTwo=fseek(fid_CardTwo,0, 'bof');  
CH3_data = fread(fid_CardTwo, inf, 'int16',skip_byte,'b');  
where_CardTwo=fseek(fid_CardTwo,2, 'bof');  
CH4_data = fread(fid_CardTwo, inf, 'int16',skip_byte,'b');  
% CH3_index=1:2:length(data_CardTwo_read);
% CH4_index=2:2:length(data_CardTwo_read);
% 
% CH3_data=data_CardTwo_read(CH3_index);
% CH4_data=data_CardTwo_read(CH4_index);



%%%触发信息和参数设置读取
data_info = importdata(strcat('H:\触发闪电\software_CardConfig_20250823172542.ini'), '', 32);
figure()
plot(CH1_data)
figure()
plot(CH2_data)
figure()
plot(CH3_data)
figure()
plot(CH4_data)


%% %%%tdms 文件读取 引雷点   
%%tdms读取方式一样，  更改文件名即可， 每个通道单独存储
string_400='H:\7天线\data\2025-08-23\20250823172542_5453CH4.tdms';  
[output_400,metaStruct_400] = TDMS_getStruct(string_400);    
int_400=output_400.groups.chans.data; %%%数据
int_trigger_t=output_400.groups.propValues;     %%%触发时间
figure
plot(int_400)
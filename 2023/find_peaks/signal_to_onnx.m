function predictedLabel = signal_to_onnx(signal, onnxPath, fs)
    % signal_to_onnx: 将输入信号转换为时间-频率图，并直接用 ONNX 模型预测类别。
    %
    % 输入参数:
    %   - signal: 输入的信号数据
    %   - onnxPath: ONNX 模型文件路径
    %   - fs: 采样频率 (Hz)
    %
    % 输出参数:
    %   - predictedLabel: ONNX 模型预测的类别标签

    % 连续小波变换 (CWT)
    [wt, f] = cwt(signal, 'amor', fs);
    time_microseconds = linspace(0, length(signal) / fs, length(signal)) * 1e6;
    frequency_MHz = f * 1e-6;

    % 转换为模型输入格式
    figure('Visible', 'off');
    imagesc(time_microseconds, frequency_MHz, abs(wt));
    axis xy;
    set(gca, 'LooseInset', get(gca, 'TightInset')); % 去除多余边框
    colormap(jet); % 使用颜色图
    frame = getframe(gca); % 获取图像帧
    close(gcf);

    % 获取图像矩阵并调整尺寸
    img = frame.cdata;
    inputImage = imresize(img, [224, 224]); % 调整为 ONNX 模型输入尺寸
    inputImage = single(inputImage) / 255; % 归一化到 [0,1]

    % 调整输入维度为 [batch_size, channels, height, width]
    inputImage = permute(inputImage, [3, 1, 2]);  % 从 [height, width, channels] 转换为 [channels, height, width]
    inputImage = reshape(inputImage, [1, size(inputImage, 1), size(inputImage, 2), size(inputImage, 3)]); % 添加 batch_size 维度

    % 使用 Python 中的 onnxruntime 推理模型
    try
        % 使用 Python onnxruntime 库加载模型
        py.importlib.import_module('onnxruntime');
        
        % 创建 InferenceSession
        session = py.onnxruntime.InferenceSession(onnxPath);
        
        % 获取输入输出节点
        inputName = session.get_inputs{1}.name;
        outputName = session.get_outputs{1}.name;
        
        % 将图像转换为 NumPy 数组格式并传入模型
        inputTensor = py.numpy.array(inputImage);
        
        % 将输入数据作为字典传递
        inputs = py.dict(pyargs(inputName, inputTensor));
        
        % 进行推理
        result = session.run({outputName}, inputs);
        [~, predictedLabel] = max(single(result{1}));  % 找到最大得分的索引，作为预测标签
        
    catch ME
        error('Error during ONNX model inference: %s', ME.message);
    end
end

clc; clear; close all;
%Parameters
M = 16;                       % QAM调制阶数
bitsPerSymbol = log2(M);      % 每个符号的比特数 (log2(16) = 4)

%修改点1: 提升N，以支持更多信道的数据传输
N = 2048;                     % OFDM点数（必须偶数，增加以支持3或4信道）
aco_carriers = 1:2:(N/2-1);   % ACO-OFDM应该用奇数子载波1,3,5,...

numTrials = 10000;            % 仿真次数（建议调试先小，后大）
snrRange = 0:0.5:14;          % 信噪比范围（dB）

% === 修改点2: 3x3 MIMO信道矩阵 ===
H = [1, 0, 0;                 % MIMO信道，对角矩阵表示独立信道
     0, 0.8, 0;
     0, 0, 0.6];              % 可以调整信道增益

% LDPC参数（WIFI常用码率，可根据实际需要替换）
% 码率0.5
[cfgEnc, cfgDec] = getProtoMatrix(648, 324); 
infoLength = cfgEnc.NumInformationBits;   % 324
blockLength = cfgEnc.BlockLength;         % 648
numIter = 25;                             % LDPC最大迭代次数

%修改点3: 信道/数据流数量 ===
streams = 3;                  

% 计算和展示实际使用的子载波数量[how many subcarriers ACO-OFDM actually uses for data, and checks the number of channels in the simulation.]
numSubcarriers = length(aco_carriers);
fprintf('使用的子载波数量: %d\n', numSubcarriers); % (N/2-1)/2 = (2048/2-1)/2 = 511.5，取整511
fprintf('当前仿真信道数量: %d\n', streams);

% ------- 结果记录 -------
% === 修改点4: 为3个数据流和总BER创建存储数组 ===
ber_interleaved_1 = zeros(1,length(snrRange));
ber_interleaved_2 = zeros(1,length(snrRange));
ber_interleaved_3 = zeros(1,length(snrRange)); 
ber_interleaved_total = zeros(1,length(snrRange));

ber_Independent_1 = zeros(1,length(snrRange));
ber_Independent_2 = zeros(1,length(snrRange));
ber_Independent_3 = zeros(1,length(snrRange)); 
ber_Independent_total = zeros(1,length(snrRange));

%  并行计算循环 
parfor snrIdx = 1:length(snrRange)
    snr = snrRange(snrIdx);
    EbN0_lin = 10^(snr/10);
    
    % 修正：考虑编码率和调制阶数的影响
    codeRate = infoLength/blockLength;  % rate= 324/648 = 0.5
    EsN0_lin = EbN0_lin * codeRate * bitsPerSymbol;  % 符号信噪比
    noiseVar = 1/(2*EsN0_lin);  % 复数噪声方差
    
    % 交织LDPC + ACO-OFDM 
    % 修改点5: error计数器扩到3路
    err1_inter = 0; err2_inter = 0; err3_inter = 0;
    
    for trial = 1:numTrials
        % 修改点6: 信息比特生成 - 均分给3路 ===
        data1 = randi([0,1], infoLength/streams, 1); % infoLength / 3
        data2 = randi([0,1], infoLength/streams, 1);
        data3 = randi([0,1], infoLength/streams, 1);
        
        % === 修改点7: 3路比特交织和解交织函数 ===
        % 交织：(a1,b1,c1,a2,b2,c2,...)
        interleave_bits = @(a,b,c) reshape([a(:)'; b(:)'; c(:)'], [], 1);
        % 解交织：分离出原始的a,b,c序列
        deinterleave_bits = @(b_total) deal(b_total(1:streams:end), b_total(2:streams:end), b_total(3:streams:end));
        
        merged = interleave_bits(data1, data2, data3); % 3路信息比特交织
        codeword = ldpcEncode(merged, cfgEnc);           % 对交织后的比特进行编码
        
        [cw1, cw2, cw3] = deinterleave_bits(codeword); % 编码后码字解交织
        
        % === 修改点8: 确保比特数能被bitsPerSymbol整除，为3路分别填充 ===
        if mod(length(cw1), bitsPerSymbol) ~= 0
            padding = bitsPerSymbol - mod(length(cw1), bitsPerSymbol);
            cw1 = [cw1; zeros(padding, 1)];
            cw2 = [cw2; zeros(padding, 1)];
            cw3 = [cw3; zeros(padding, 1)]; % 第三路填充
        end
        
        % === 修改点9: QAM调制3路比特流 ===
        symbols1 = qammod(cw1, M, 'InputType','bit','UnitAveragePower',true);
        symbols2 = qammod(cw2, M, 'InputType','bit','UnitAveragePower',true);
        symbols3 = qammod(cw3, M, 'InputType','bit','UnitAveragePower',true); % 第三路QAM调制
        
        % 确保符号数不超过可用子载波数
        % (现在N=2048，numSubcarriers=511，每路QAM符号数为 blockLength/bitsPerSymbol = 648/4 = 162，
        %  合并后总符号数162，但交织后符号数是 blockLength/bitsPerSymbol = 648/4 = 162，
        %  这个值远小于 511，所以不会发生截断)
        numSymbolsPerStream = min(length(symbols1), numSubcarriers);
        symbols1 = symbols1(1:numSymbolsPerStream);
        symbols2 = symbols2(1:numSymbolsPerStream);
        symbols3 = symbols3(1:numSymbolsPerStream); % 第三路符号截取
        aco_carriers_use = aco_carriers(1:numSymbolsPerStream);
        
        % === 修改点10: ACO-OFDM映射到3路频域信号 ===
        X1 = zeros(N,1); X2 = zeros(N,1); X3 = zeros(N,1); % 频域OFDM符号初始化
        X1(aco_carriers_use+1) = symbols1;  % 映射符号到奇数子载波
        X2(aco_carriers_use+1) = symbols2;
        X3(aco_carriers_use+1) = symbols3;
        
        % Hermitian对称性确保IFFT后实信号
        X1(N-aco_carriers_use+1) = conj(symbols1);
        X2(N-aco_carriers_use+1) = conj(symbols2);
        X3(N-aco_carriers_use+1) = conj(symbols3);
        
        % === 修改点11: IFFT并添加ACO剪切（3路） ===
        x1 = ifft(X1)*sqrt(N); x2 = ifft(X2)*sqrt(N); x3 = ifft(X3)*sqrt(N);
        % ACO：负值置零，正值乘以2保持功率
        x1_aco = 2*max(real(x1), 0);
        x2_aco = 2*max(real(x2), 0);
        x3_aco = 2*max(real(x3), 0);
        
        % === 修改点12: MIMO信道传输 ===
        tx = [x1_aco x2_aco x3_aco]; % 发送信号矩阵 (N x 3)
        noise = sqrt(noiseVar)*(randn(N,streams) + 1j*randn(N,streams))/sqrt(2); % 噪声维度匹配 (N x 3)
        rx = tx*H + noise; % 接收信号[awgn]
        
        % === 修改点13: FFT/子载波提取/等化（3路） ===
        Y1 = fft(rx(:,1))/sqrt(N); Y2 = fft(rx(:,2))/sqrt(N); Y3 = fft(rx(:,3))/sqrt(N);
        rxsyms1 = Y1(aco_carriers_use+1)/H(1,1); % 等化第一路
        rxsyms2 = Y2(aco_carriers_use+1)/H(2,2); % 等化第二路
        rxsyms3 = Y3(aco_carriers_use+1)/H(3,3); % 等化第三路
        
        % === 修改点14: 软解调并截取到原始长度（3路） ===
        effectiveNoiseVar1 = noiseVar/(H(1,1)^2);
        effectiveNoiseVar2 = noiseVar/(H(2,2)^2);
        effectiveNoiseVar3 = noiseVar/(H(3,3)^2); % 第三路有效噪声方差
        
        llr1 = qamdemod(rxsyms1, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar1);
        llr2 = qamdemod(rxsyms2, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar2);
        llr3 = qamdemod(rxsyms3, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar3); % 第三路软解调
        
        llr1 = llr1(1:length(cw1)); % 截取到原始编码比特长度
        llr2 = llr2(1:length(cw2));
        llr3 = llr3(1:length(cw3));
        
        % === 修改点15: LLR交织和LDPC解码 ===
        llr_total = interleave_bits(llr1, llr2, llr3); % LLR进行交织
        decoded_total = ldpcDecode(llr_total, cfgDec, numIter); % 对交织后的LLR进行解码
        
        [dec1, dec2, dec3] = deinterleave_bits(decoded_total); % 解码结果解交织
        
        % === 修改点16: 错误比特计数（只比较原始信息位） ===
        err1_inter = err1_inter + sum(dec1(1:length(data1)) ~= data1);
        err2_inter = err2_inter + sum(dec2(1:length(data2)) ~= data2);
        err3_inter = err3_inter + sum(dec3(1:length(data3)) ~= data3); % 第三路错误计数
    end
    
    % === 修改点17: 计算误码率并存储（交织LDPC） ===
    ber_interleaved_1(snrIdx) = err1_inter / (numTrials * infoLength/streams);
    ber_interleaved_2(snrIdx) = err2_inter / (numTrials * infoLength/streams);
    ber_interleaved_3(snrIdx) = err3_inter / (numTrials * infoLength/streams); % 第三路BER
    ber_interleaved_total(snrIdx) = (err1_inter + err2_inter + err3_inter) / (numTrials * infoLength);
    
    %% 独立LDPC + ACO-OFDM
    %修改点18: error计数器扩展到3路 ===
    err1_Independent = 0; err2_Independent = 0; err3_Independent = 0;
    
    for trial = 1:numTrials
        %  修改点19: 信息比特生成（每路独立编码，所以每个是infoLength） ===
        bits1 = randi([0,1], infoLength, 1);
        bits2 = randi([0,1], infoLength, 1);
        bits3 = randi([0,1], infoLength, 1); % 第三路原始比特
        
        % === 修改点20: 独立LDPC编码（3路） ===
        cw1 = ldpcEncode(bits1, cfgEnc);
        cw2 = ldpcEncode(bits2, cfgEnc);
        cw3 = ldpcEncode(bits3, cfgEnc); % 第三路独立编码
        
        % === 修改点21: 确保比特数能被bitsPerSymbol整除，为3路分别填充 ===
        if mod(length(cw1), bitsPerSymbol) ~= 0
            padding = bitsPerSymbol - mod(length(cw1), bitsPerSymbol);
            cw1 = [cw1; zeros(padding, 1)];
            cw2 = [cw2; zeros(padding, 1)];
            cw3 = [cw3; zeros(padding, 1)]; % 第三路填充
        end
        
        % === 修改点22: QAM调制3路比特流 ===
        symbols1 = qammod(cw1, M, 'InputType','bit','UnitAveragePower',true);
        symbols2 = qammod(cw2, M, 'InputType','bit','UnitAveragePower',true);
        symbols3 = qammod(cw3, M, 'InputType','bit','UnitAveragePower',true); % 第三路QAM调制
        
        % 确保符号数不超过可用子载波数
        % (现在N=2048，numSubcarriers=511，每路QAM符号数为 blockLength/bitsPerSymbol = 648/4 = 162，
        %  总符号数 162*3 = 486，小于511，所以不会发生截断)
        numSymbolsPerStream = min(length(symbols1), numSubcarriers);
        symbols1 = symbols1(1:numSymbolsPerStream);
        symbols2 = symbols2(1:numSymbolsPerStream);
        symbols3 = symbols3(1:numSymbolsPerStream);
        aco_carriers_use = aco_carriers(1:numSymbolsPerStream);
        
        % === 修改点23: ACO-OFDM映射到3路频域信号 ===
        X1 = zeros(N,1); X2 = zeros(N,1); X3 = zeros(N,1);
        X1(aco_carriers_use+1) = symbols1;
        X2(aco_carriers_use+1) = symbols2;
        X3(aco_carriers_use+1) = symbols3;
        X1(N-aco_carriers_use+1) = conj(symbols1);
        X2(N-aco_carriers_use+1) = conj(symbols2);
        X3(N-aco_carriers_use+1) = conj(symbols3);
        
        % === 修改点24: IFFT并添加ACO剪切（3路） ===
        x1 = ifft(X1)*sqrt(N); x2 = ifft(X2)*sqrt(N); x3 = ifft(X3)*sqrt(N);
        x1_aco = 2*max(real(x1), 0);
        x2_aco = 2*max(real(x2), 0);
        x3_aco = 2*max(real(x3), 0);
        
        % === 修改点25: MIMO信道传输 ===
        tx = [x1_aco x2_aco x3_aco];
        noise = sqrt(noiseVar)*(randn(N,streams) + 1j*randn(N,streams))/sqrt(2);
        rx = tx*H + noise;
        
        % === 修改点26: FFT/子载波提取/等化（3路） ===
        Y1 = fft(rx(:,1))/sqrt(N); Y2 = fft(rx(:,2))/sqrt(N); Y3 = fft(rx(:,3))/sqrt(N);
        rxsyms1 = Y1(aco_carriers_use+1)/H(1,1);
        rxsyms2 = Y2(aco_carriers_use+1)/H(2,2);
        rxsyms3 = Y3(aco_carriers_use+1)/H(3,3);
        
        % === 修改点27: 软解调并截取到原始长度（3路） ===
        effectiveNoiseVar1 = noiseVar/(H(1,1)^2);
        effectiveNoiseVar2 = noiseVar/(H(2,2)^2);
        effectiveNoiseVar3 = noiseVar/(H(3,3)^2);
        
        llr1 = qamdemod(rxsyms1, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar1);
        llr2 = qamdemod(rxsyms2, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar2);
        llr3 = qamdemod(rxsyms3, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar3);
        
        llr1 = llr1(1:length(cw1));
        llr2 = llr2(1:length(cw2));
        llr3 = llr3(1:length(cw3));
        
        % === 修改点28: 独立LDPC解码（3路） ===
        dec1 = ldpcDecode(llr1, cfgDec, numIter);
        dec2 = ldpcDecode(llr2, cfgDec, numIter);
        dec3 = ldpcDecode(llr3, cfgDec, numIter); % 第三路独立解码
        
        % === 修改点29: 错误比特计数（只比较原始信息位） ===
        err1_Independent = err1_Independent + sum(dec1(1:length(bits1)) ~= bits1);
        err2_Independent = err2_Independent + sum(dec2(1:length(bits2)) ~= bits2);
        err3_Independent = err3_Independent + sum(dec3(1:length(bits3)) ~= bits3); % 第三路错误计数
    end
    
    % === 修改点30: 计算误码率并存储（独立LDPC） ===
    ber_Independent_1(snrIdx) = err1_Independent / (numTrials * infoLength);
    ber_Independent_2(snrIdx) = err2_Independent / (numTrials * infoLength);
    ber_Independent_3(snrIdx) = err3_Independent / (numTrials * infoLength); % 第三路BER
    ber_Independent_total(snrIdx) = (err1_Independent + err2_Independent + err3_Independent) / (numTrials * infoLength * streams); % 总发送比特数 = numTrials * infoLength * 3
    
    % === 修改点31: 命令行打印结果 ===
    fprintf('SNR=%2d: 交织LDPC S1=%.2e S2=%.2e S3=%.2e | 独立LDPC S1=%.2e S2=%.2e S3=%.2e\n', ...
        snr, ber_interleaved_1(snrIdx), ber_interleaved_2(snrIdx), ber_interleaved_3(snrIdx), ...
        ber_Independent_1(snrIdx), ber_Independent_2(snrIdx), ber_Independent_3(snrIdx));
end

% ------- 画图 -------
figure('Position', [100, 100, 1200, 600]);

% 交织LDPC BER曲线 ===
subplot(1,2,1);
semilogy(snrRange, ber_Independent_total, 'r-s', 'LineWidth', 2, 'DisplayName','Overall BER');
hold on;
semilogy(snrRange, ber_Independent_1, 'y-d', 'LineWidth', 1.5, 'DisplayName','Stream 1 (H=1)');
semilogy(snrRange, ber_Independent_2, 'k-p', 'LineWidth', 1.5, 'DisplayName','Stream 2 (H=0.8)');
semilogy(snrRange, ber_Independent_3, 'b-h', 'LineWidth', 1.5, 'DisplayName','Stream 3 (H=0.6)'); % 新增第三路曲线
grid on;
xlabel('Eb/N0 (dB)','FontSize', 10); 
ylabel('Bit Error Rate (BER)','FontSize', 14);
title({'Independent LDPC (3-Stream)'},'FontSize', 14);
legend('show', 'Location', 'southwest','FontSize', 14);
ylim([1e-6 1]);xlim([0 6]); xticks(0:0.5:6);
ax1 = gca; ax1.FontSize = 14; 
%  独立LDPC BER曲线 ===
subplot(1,2,2);
semilogy(snrRange, ber_interleaved_total, 'b-o', 'LineWidth', 2, 'DisplayName','Overall BER');
hold on;
semilogy(snrRange, ber_interleaved_1, 'c--^', 'LineWidth', 1.5, 'DisplayName','Stream 1 (H=1)');
semilogy(snrRange, ber_interleaved_2, 'm--v', 'LineWidth', 1.5, 'DisplayName','Stream 2 (H=0.8)');
semilogy(snrRange, ber_interleaved_3, 'g--x', 'LineWidth', 1.5, 'DisplayName','Stream 3 (H=0.6)'); % 新增第三路曲线
grid on;
xlabel('Eb/N0 (dB)','FontSize', 10); 
ylabel('Bit Error Rate (BER)','FontSize', 14);
title({'interleaved LDPC (3-Stream)'},'FontSize', 14);
legend('show', 'Location', 'southwest','FontSize', 14); 
ylim([1e-6 1]);
xlim([0 6]); xticks(0:0.5:6);  
ax2 = gca; 
ax2.FontSize = 14; 
%compare
figure('Position', [100, 100, 400, 500]); 
semilogy(snrRange, ber_Independent_total, 'r-s', 'LineWidth', 2, 'DisplayName','Independent LDPC Overall BER');
hold on;
semilogy(snrRange, ber_interleaved_total, 'b-o', 'LineWidth', 2, 'DisplayName','Interleaved LDPC Overall BER');
grid on;
xlabel('Eb/N0 (dB)', 'FontSize', 10);
ylabel('Bit Error Rate (BER)', 'FontSize', 14);
title({'Overall BER Performance'}, 'FontSize', 14); 
legend('show', 'Location', 'southwest', 'FontSize', 14);
ylim([1e-6 1]); 
xlim([0 6]);   
xticks(0:0.5:6); 
ax3 = gca; 
ax3.FontSize = 14; 
clc; clear; close all;
fprintf('--- MIMO ACO-OFDM: Interleaved LDPC vs Dual-LDPC 对比 ---\n');

% 基本参数
M = 16;                       % QAM调制阶数
bitsPerSymbol = log2(M);
N = 1024;                     % OFDM点数（必须偶数）
aco_carriers = 1:2:(N/2-1);   % 修正：ACO-OFDM应该用奇数子载波1,3,5,...
numTrials = 10000;              % 仿真次数（建议调试先小，后大）
snrRange = 0:0.5:14;            % 信噪比范围（dB）
H = [1, 0; 
     0, 0.8];           % MIMO信道

% LDPC参数（这里是WIFI常用码率，可根据实际需要替换）
[cfgEnc, cfgDec] = getProtoMatrix(648, 324); % 648码长, 324信息位
infoLength = cfgEnc.NumInformationBits;       % 324
blockLength = cfgEnc.BlockLength;             % 648

numIter = 25;                  % LDPC最大迭代
streams = 2;                   % 双路

% 计算实际使用的子载波数量
numSubcarriers = length(aco_carriers);
fprintf('使用的子载波数量: %d\n', numSubcarriers);

% ------- 结果记录 -------
ber_interleaved_1 = zeros(1,length(snrRange));
ber_interleaved_2 = zeros(1,length(snrRange));
ber_interleaved_total = zeros(1,length(snrRange));

ber_dual_1 = zeros(1,length(snrRange));
ber_dual_2 = zeros(1,length(snrRange));
ber_dual_total = zeros(1,length(snrRange));

parfor snrIdx = 1:length(snrRange)
    snr = snrRange(snrIdx);
    EbN0_lin = 10^(snr/10);
    % 修正：考虑编码率和调制阶数的影响
    codeRate = infoLength/blockLength;  % 编码率 = 324/648 = 0.5
    EsN0_lin = EbN0_lin * codeRate * bitsPerSymbol;  % 符号信噪比
    noiseVar = 1/(2*EsN0_lin);  % 复数噪声方差

    %交织LDPC + ACO-OFDM 
    err1 = 0; err2 = 0;
    for trial = 1:numTrials
        % 信息比特
        data1 = randi([0,1], infoLength/2, 1);
        data2 = randi([0,1], infoLength/2, 1);
        % 交织、编码
        interleave_bits = @(a,b) reshape([a(:)'; b(:)'], [], 1);
        deinterleave_bits = @(b) deal(b(1:2:end), b(2:2:end));
        merged = interleave_bits(data1, data2);
        codeword = ldpcEncode(merged, cfgEnc);
        [cw1, cw2] = deinterleave_bits(codeword);
        
        % 修正：确保比特数能被bitsPerSymbol整除
        if mod(length(cw1), bitsPerSymbol) ~= 0
            padding = bitsPerSymbol - mod(length(cw1), bitsPerSymbol);
            cw1 = [cw1; zeros(padding, 1)];
            cw2 = [cw2; zeros(padding, 1)];
        end
        
        % QAM调制
        symbols1 = qammod(cw1, M, 'InputType','bit','UnitAveragePower',true);
        symbols2 = qammod(cw2, M, 'InputType','bit','UnitAveragePower',true);
        
        % 确保符号数不超过可用子载波数
        numSymbolsPerStream = min(length(symbols1), numSubcarriers);
        symbols1 = symbols1(1:numSymbolsPerStream);
        symbols2 = symbols2(1:numSymbolsPerStream);
        aco_carriers_use = aco_carriers(1:numSymbolsPerStream);
        
        % 填子载波 - ACO-OFDM的正确映射
        X1 = zeros(N,1); X2 = zeros(N,1);
        X1(aco_carriers_use+1) = symbols1;  % 修正：MATLAB索引从1开始
        X2(aco_carriers_use+1) = symbols2;
        % Hermitian对称性确保实信号
        X1(N-aco_carriers_use+1) = conj(symbols1);
        X2(N-aco_carriers_use+1) = conj(symbols2);
        
        % IFFT并添加ACO剪切
        x1 = ifft(X1)*sqrt(N); x2 = ifft(X2)*sqrt(N);
        % ACO剪切：负值置零，正值乘以2保持功率
        x1_aco = 2*max(real(x1), 0); x2_aco = 2*max(real(x2), 0);
        
        tx = [x1_aco x2_aco];
        noise = sqrt(noiseVar)*(randn(N,2) + 1j*randn(N,2))/sqrt(2);
        rx = tx*H + noise;
        
        % FFT/子载波提取/等化
        Y1 = fft(rx(:,1))/sqrt(N); Y2 = fft(rx(:,2))/sqrt(N);
        rxsyms1 = Y1(aco_carriers_use+1)/H(1,1);
        rxsyms2 = Y2(aco_carriers_use+1)/H(2,2);
        
        % 软解调 - 修正噪声方差
        effectiveNoiseVar1 = noiseVar/(H(1,1)^2);
        effectiveNoiseVar2 = noiseVar/(H(2,2)^2);
        
        llr1 = qamdemod(rxsyms1, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar1);
        llr2 = qamdemod(rxsyms2, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar2);
        
        % 截取到原始长度
        llr1 = llr1(1:length(cw1));
        llr2 = llr2(1:length(cw2));
        
        llr_total = interleave_bits(llr1, llr2);
        decoded = ldpcDecode(llr_total, cfgDec, numIter);
        [dec1, dec2] = deinterleave_bits(decoded);
        
        % 只比较原始数据位
        err1 = err1 + sum(dec1(1:length(data1)) ~= data1);
        err2 = err2 + sum(dec2(1:length(data2)) ~= data2);
    end
    ber_interleaved_1(snrIdx) = err1 / (numTrials * infoLength/2);
    ber_interleaved_2(snrIdx) = err2 / (numTrials * infoLength/2);
    ber_interleaved_total(snrIdx) = (err1 + err2) / (numTrials * infoLength);

    %双路LDPC + ACO-OFDM -----------
    err1 = 0; err2 = 0;
    for trial = 1:numTrials
        bits1 = randi([0,1], infoLength, 1);
        bits2 = randi([0,1], infoLength, 1);
        cw1 = ldpcEncode(bits1, cfgEnc);
        cw2 = ldpcEncode(bits2, cfgEnc);
        
        % 修正：确保比特数能被bitsPerSymbol整除
        if mod(length(cw1), bitsPerSymbol) ~= 0
            padding = bitsPerSymbol - mod(length(cw1), bitsPerSymbol);
            cw1 = [cw1; zeros(padding, 1)];
            cw2 = [cw2; zeros(padding, 1)];
        end
        
        symbols1 = qammod(cw1, M, 'InputType','bit','UnitAveragePower',true);
        symbols2 = qammod(cw2, M, 'InputType','bit','UnitAveragePower',true);
        
        % 确保符号数不超过可用子载波数
        numSymbolsPerStream = min(length(symbols1), numSubcarriers);
        symbols1 = symbols1(1:numSymbolsPerStream);
        symbols2 = symbols2(1:numSymbolsPerStream);
        aco_carriers_use = aco_carriers(1:numSymbolsPerStream);
        
        % ACO-OFDM映射
        X1 = zeros(N,1); X2 = zeros(N,1);
        X1(aco_carriers_use+1) = symbols1;
        X2(aco_carriers_use+1) = symbols2;
        X1(N-aco_carriers_use+1) = conj(symbols1);
        X2(N-aco_carriers_use+1) = conj(symbols2);
        
        x1 = ifft(X1)*sqrt(N); x2 = ifft(X2)*sqrt(N);
        x1_aco = 2*max(real(x1), 0); x2_aco = 2*max(real(x2), 0);
        
        tx = [x1_aco x2_aco];
        noise = sqrt(noiseVar)*(randn(N,2) + 1j*randn(N,2))/sqrt(2);
        rx = tx*H + noise;
        
        Y1 = fft(rx(:,1))/sqrt(N); Y2 = fft(rx(:,2))/sqrt(N);
        rxsyms1 = Y1(aco_carriers_use+1)/H(1,1);
        rxsyms2 = Y2(aco_carriers_use+1)/H(2,2);
        
        effectiveNoiseVar1 = noiseVar/(H(1,1)^2);
        effectiveNoiseVar2 = noiseVar/(H(2,2)^2);
        
        llr1 = qamdemod(rxsyms1, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar1);
        llr2 = qamdemod(rxsyms2, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar2);
        
        % 截取到原始长度
        llr1 = llr1(1:length(cw1));
        llr2 = llr2(1:length(cw2));
        
        dec1 = ldpcDecode(llr1, cfgDec, numIter);
        dec2 = ldpcDecode(llr2, cfgDec, numIter);
        
        err1 = err1 + sum(dec1(1:length(bits1)) ~= bits1);
        err2 = err2 + sum(dec2(1:length(bits2)) ~= bits2);
    end
    ber_dual_1(snrIdx) = err1 / (numTrials * infoLength);
    ber_dual_2(snrIdx) = err2 / (numTrials * infoLength);
    ber_dual_total(snrIdx) = (err1 + err2) / (numTrials * infoLength * 2);

    fprintf('SNR=%2d: 交织LDPC S1=%.2e S2=%.2e | 双LDPC S1=%.2e S2=%.2e\n', ...
        snr, ber_interleaved_1(snrIdx), ber_interleaved_2(snrIdx), ...
        ber_dual_1(snrIdx), ber_dual_2(snrIdx));
end

% 画图 
figure('Position', [100, 100, 1200, 600]);
subplot(1,2,1);
semilogy(snrRange, ber_interleaved_total, 'b-o', 'LineWidth', 2, 'DisplayName','Overall BER');
hold on;
semilogy(snrRange, ber_interleaved_1, 'c--^', 'LineWidth', 1.5, 'DisplayName','Stream 1 (H=1)');
semilogy(snrRange, ber_interleaved_2, 'm--v', 'LineWidth', 1.5, 'DisplayName','Stream 2 (H=0.8)');
grid on;
xlabel('Eb/N0 (dB)'); ylabel('Bit Error Rate (BER)');
title({'Architecture 1: Interleaved LDPC'});
legend('show', 'Location', 'southwest'); ylim([1e-6 1]);

subplot(1,2,2);
semilogy(snrRange, ber_dual_total, 'r-s', 'LineWidth', 2, 'DisplayName','Overall BER');
hold on;
semilogy(snrRange, ber_dual_1, 'y-d', 'LineWidth', 1.5, 'DisplayName','Stream 1 (H=1)');
semilogy(snrRange, ber_dual_2, 'k-p', 'LineWidth', 1.5, 'DisplayName','Stream 2 (H=0.8)');
grid on;
xlabel('Eb/N0 (dB)'); ylabel('Bit Error Rate (BER)');
title({'Architecture 2: Original Dual LDPC'});
legend('show', 'Location', 'southwest'); ylim([1e-6 1]);

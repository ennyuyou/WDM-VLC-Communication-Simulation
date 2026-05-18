clc; clear; close all;
fprintf('--- MIMO ACO-OFDM: 4信道交织LDPC vs 独立LDPC 对比 (自定义信道增益) ---\n');

% 参数
M = 16;                       % 16-QAM
bitsPerSymbol = log2(M);      %log2(16) = 4

N = 4096;                     % OFDM点数 (FFT/IFFT长度)。必须是偶数。
                              % 增大N是为了确保有足够的奇数子载波支持多路传输。
aco_carriers = 1:2:(N/2-1);   % ACO-OFDM仅使用奇数子载波进行数据传输 (1, 3, 5, ..., N/2-1)

numTrials = 10000;            % 仿真次数：每个信噪比点重复发送和接收数据的次数，用于统计误码率。
snrRange = 0:0.5:14;          % 信噪比(Eb/N0)范围（单位：dB）。

streams = 4;                  % 系统中的数据流/MIMO信道数量。

H = [1, 0, 0, 0;              
     0, 0.9, 0, 0;  
     0, 0, 0.5, 0;  
     0, 0, 0, 0.2]; 

% LDPC编码参数
[cfgEnc, cfgDec] = getProtoMatrix(648, 324); 
infoLength = cfgEnc.NumInformationBits;   
blockLength = cfgEnc.BlockLength;         
numIter = 25;                             

if isprop(cfgDec, 'MaxNumIterations')
    cfgDec.MaxNumIterations = numIter; 
end

numSubcarriers = length(aco_carriers);

fprintf('使用的OFDM点数 (N): %d\n', N);
fprintf('可用的奇数子载波数量: %d\n', numSubcarriers); 
fprintf('当前仿真MIMO信道数量: %d\n', streams);
fprintf('LDPC码字长: %d, 信息位长: %d, 码率: %.1f\n', blockLength, infoLength, infoLength/blockLength);
fprintf('自定义信道增益 (H对角线): [1, 0.9, 0.5, 0.2]\n');


%结果记录 
ber_interleaved_1 = zeros(1,length(snrRange));
ber_interleaved_2 = zeros(1,length(snrRange));
ber_interleaved_3 = zeros(1,length(snrRange));
ber_interleaved_4 = zeros(1,length(snrRange));
ber_interleaved_total = zeros(1,length(snrRange));

ber_dual_1 = zeros(1,length(snrRange));
ber_dual_2 = zeros(1,length(snrRange));
ber_dual_3 = zeros(1,length(snrRange));
ber_dual_4 = zeros(1,length(snrRange));
ber_dual_total = zeros(1,length(snrRange));

% 主仿真循环：遍历不同的信噪比 ===
parfor snrIdx = 1:length(snrRange)
    snr = snrRange(snrIdx);
    EbN0_lin = 10^(snr/10); 
    
    codeRate = infoLength/blockLength;  
    EsN0_lin = EbN0_lin * codeRate * bitsPerSymbol;  
    noiseVar = 1/(2*EsN0_lin);  
    
    %%架构1: 交织LDPC + MIMO ACO-OFDM 仿真
    err1_inter = 0; err2_inter = 0; err3_inter = 0; err4_inter = 0;
    
    for trial = 1:numTrials
        data1 = randi([0,1], infoLength/streams, 1); 
        data2 = randi([0,1], infoLength/streams, 1);
        data3 = randi([0,1], infoLength/streams, 1);
        data4 = randi([0,1], infoLength/streams, 1);
        
        interleave_bits = @(a,b,c,d) reshape([a(:)'; b(:)'; c(:)'; d(:)'], [], 1);
        merged_info = interleave_bits(data1, data2, data3, data4); 

        codeword = ldpcEncode(merged_info, cfgEnc);                  
        
        deinterleave_bits = @(b_total) deal(b_total(1:streams:end), b_total(2:streams:end), b_total(3:streams:end), b_total(4:streams:end));
        [cw1, cw2, cw3, cw4] = deinterleave_bits(codeword); 
        
        if mod(length(cw1), bitsPerSymbol) ~= 0
            padding = bitsPerSymbol - mod(length(cw1), bitsPerSymbol);
            cw1 = [cw1; zeros(padding, 1)];
            cw2 = [cw2; zeros(padding, 1)];
            cw3 = [cw3; zeros(padding, 1)];
            cw4 = [cw4; zeros(padding, 1)]; 
        end
        
        symbols1 = qammod(cw1, M, 'InputType','bit','UnitAveragePower',true);
        symbols2 = qammod(cw2, M, 'InputType','bit','UnitAveragePower',true);
        symbols3 = qammod(cw3, M, 'InputType','bit','UnitAveragePower',true);
        symbols4 = qammod(cw4, M, 'InputType','bit','UnitAveragePower',true); 
        
        numSymbolsPerStream = min(length(symbols1), numSubcarriers);
        symbols1 = symbols1(1:numSymbolsPerStream);
        symbols2 = symbols2(1:numSymbolsPerStream);
        symbols3 = symbols3(1:numSymbolsPerStream);
        symbols4 = symbols4(1:numSymbolsPerStream);
        aco_carriers_use = aco_carriers(1:numSymbolsPerStream);

        X1 = zeros(N,1); X2 = zeros(N,1); X3 = zeros(N,1); X4 = zeros(N,1); 
        X1(aco_carriers_use+1) = symbols1;  
        X2(aco_carriers_use+1) = symbols2;
        X3(aco_carriers_use+1) = symbols3;
        X4(aco_carriers_use+1) = symbols4;
        
        X1(N-aco_carriers_use+1) = conj(symbols1);
        X2(N-aco_carriers_use+1) = conj(symbols2);
        X3(N-aco_carriers_use+1) = conj(symbols3);
        X4(N-aco_carriers_use+1) = conj(symbols4);
        
        x1 = ifft(X1)*sqrt(N); x2 = ifft(X2)*sqrt(N); x3 = ifft(X3)*sqrt(N); x4 = ifft(X4)*sqrt(N);
        x1_aco = 2*max(real(x1), 0);
        x2_aco = 2*max(real(x2), 0);
        x3_aco = 2*max(real(x3), 0);
        x4_aco = 2*max(real(x4), 0);
        
        tx = [x1_aco x2_aco x3_aco x4_aco]; 
        noise = sqrt(noiseVar)*(randn(N,streams) + 1j*randn(N,streams))/sqrt(2); 
        rx = tx*H + noise; 
        
        Y1 = fft(rx(:,1))/sqrt(N); Y2 = fft(rx(:,2))/sqrt(N); Y3 = fft(rx(:,3))/sqrt(N); Y4 = fft(rx(:,4))/sqrt(N);
        rxsyms1 = Y1(aco_carriers_use+1)/H(1,1); 
        rxsyms2 = Y2(aco_carriers_use+1)/H(2,2); 
        rxsyms3 = Y3(aco_carriers_use+1)/H(3,3); 
        rxsyms4 = Y4(aco_carriers_use+1)/H(4,4); 
        
        effectiveNoiseVar1 = noiseVar/(H(1,1)^2);
        effectiveNoiseVar2 = noiseVar/(H(2,2)^2);
        effectiveNoiseVar3 = noiseVar/(H(3,3)^2);
        effectiveNoiseVar4 = noiseVar/(H(4,4)^2); 
        
        llr1 = qamdemod(rxsyms1, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar1);
        llr2 = qamdemod(rxsyms2, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar2);
        llr3 = qamdemod(rxsyms3, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar3);
        llr4 = qamdemod(rxsyms4, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar4);
        
        llr1 = llr1(1:blockLength/streams); 
        llr1 = double(llr1(:)); 
        llr2 = llr2(1:blockLength/streams); 
        llr2 = double(llr2(:)); 
        llr3 = llr3(1:blockLength/streams); 
        llr3 = double(llr3(:)); 
        llr4 = llr4(1:blockLength/streams); 
        llr4 = double(llr4(:)); 
        
        llr_total = interleave_bits(llr1, llr2, llr3, llr4); 
        llr_total = double(llr_total(:)); 
        
        [decoded_total,~,~] = ldpcDecode(llr_total, cfgDec, numIter); 
        
        [dec1, dec2, dec3, dec4] = deinterleave_bits(decoded_total); 
        
        err1_inter = err1_inter + sum(dec1(1:length(data1)) ~= data1);
        err2_inter = err2_inter + sum(dec2(1:length(data2)) ~= data2);
        err3_inter = err3_inter + sum(dec3(1:length(data3)) ~= data3);
        err4_inter = err4_inter + sum(dec4(1:length(data4)) ~= data4); 
    end
    
    ber_interleaved_1(snrIdx) = err1_inter / (numTrials * infoLength/streams);
    ber_interleaved_2(snrIdx) = err2_inter / (numTrials * infoLength/streams);
    ber_interleaved_3(snrIdx) = err3_inter / (numTrials * infoLength/streams);
    ber_interleaved_4(snrIdx) = err4_inter / (numTrials * infoLength/streams); 
    ber_interleaved_total(snrIdx) = (err1_inter + err2_inter + err3_inter + err4_inter) / (numTrials * infoLength);
    
    %%架构2: 独立LDPC + MIMO ACO-OFDM 仿真块 
    err1_dual = 0; err2_dual = 0; err3_dual = 0; err4_dual = 0;
    
    for trial = 1:numTrials
        bits1 = randi([0,1], infoLength, 1);
        bits2 = randi([0,1], infoLength, 1);
        bits3 = randi([0,1], infoLength, 1);
        bits4 = randi([0,1], infoLength, 1); 
        
        cw1 = ldpcEncode(bits1, cfgEnc);
        cw2 = ldpcEncode(bits2, cfgEnc);
        cw3 = ldpcEncode(bits3, cfgEnc);
        cw4 = ldpcEncode(bits4, cfgEnc); 
        
        if mod(length(cw1), bitsPerSymbol) ~= 0
            padding = bitsPerSymbol - mod(length(cw1), bitsPerSymbol);
            cw1 = [cw1; zeros(padding, 1)];
            cw2 = [cw2; zeros(padding, 1)];
            cw3 = [cw3; zeros(padding, 1)];
            cw4 = [cw4; zeros(padding, 1)]; 
        end
        
        symbols1 = qammod(cw1, M, 'InputType','bit','UnitAveragePower',true);
        symbols2 = qammod(cw2, M, 'InputType','bit','UnitAveragePower',true);
        symbols3 = qammod(cw3, M, 'InputType','bit','UnitAveragePower',true);
        symbols4 = qammod(cw4, M, 'InputType','bit','UnitAveragePower',true); 
        
        numSymbolsPerStream = min(length(symbols1), numSubcarriers);
        symbols1 = symbols1(1:numSymbolsPerStream);
        symbols2 = symbols2(1:numSymbolsPerStream);
        symbols3 = symbols3(1:numSymbolsPerStream);
        symbols4 = symbols4(1:numSymbolsPerStream);
        aco_carriers_use = aco_carriers(1:numSymbolsPerStream);
        
        X1 = zeros(N,1); X2 = zeros(N,1); X3 = zeros(N,1); X4 = zeros(N,1);
        X1(aco_carriers_use+1) = symbols1;
        X2(aco_carriers_use+1) = symbols2;
        X3(aco_carriers_use+1) = symbols3;
        X4(aco_carriers_use+1) = symbols4;
        X1(N-aco_carriers_use+1) = conj(symbols1);
        X2(N-aco_carriers_use+1) = conj(symbols2);
        X3(N-aco_carriers_use+1) = conj(symbols3);
        X4(N-aco_carriers_use+1) = conj(symbols4);
        
        x1 = ifft(X1)*sqrt(N); x2 = ifft(X2)*sqrt(N); x3 = ifft(X3)*sqrt(N); x4 = ifft(X4)*sqrt(N);
        x1_aco = 2*max(real(x1), 0);
        x2_aco = 2*max(real(x2), 0);
        x3_aco = 2*max(real(x3), 0);
        x4_aco = 2*max(real(x4), 0);
        
        tx = [x1_aco x2_aco x3_aco x4_aco];
        noise = sqrt(noiseVar)*(randn(N,streams) + 1j*randn(N,streams))/sqrt(2);
        rx = tx*H + noise;
        
        Y1 = fft(rx(:,1))/sqrt(N); Y2 = fft(rx(:,2))/sqrt(N); Y3 = fft(rx(:,3))/sqrt(N); Y4 = fft(rx(:,4))/sqrt(N);
        rxsyms1 = Y1(aco_carriers_use+1)/H(1,1);
        rxsyms2 = Y2(aco_carriers_use+1)/H(2,2);
        rxsyms3 = Y3(aco_carriers_use+1)/H(3,3);
        rxsyms4 = Y4(aco_carriers_use+1)/H(4,4);
        
        effectiveNoiseVar1 = noiseVar/(H(1,1)^2);
        effectiveNoiseVar2 = noiseVar/(H(2,2)^2);
        effectiveNoiseVar3 = noiseVar/(H(3,3)^2);
        effectiveNoiseVar4 = noiseVar/(H(4,4)^2);
        
        llr1 = qamdemod(rxsyms1, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar1);
        llr2 = qamdemod(rxsyms2, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar2);
        llr3 = qamdemod(rxsyms3, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar3);
        llr4 = qamdemod(rxsyms4, M, 'OutputType','approxllr','UnitAveragePower',true, 'NoiseVariance', effectiveNoiseVar4);
        
        llr1 = llr1(1:blockLength); 
        llr1 = double(llr1(:)); 
        llr2 = llr2(1:blockLength); 
        llr2 = double(llr2(:)); 
        llr3 = llr3(1:blockLength); 
        llr3 = double(llr3(:)); 
        llr4 = llr4(1:blockLength); 
        llr4 = double(llr4(:)); 
        
        [decoded1,~,~] = ldpcDecode(llr1, cfgDec, numIter);
        [decoded2,~,~] = ldpcDecode(llr2, cfgDec, numIter);
        [decoded3,~,~] = ldpcDecode(llr3, cfgDec, numIter);
        [decoded4,~,~] = ldpcDecode(llr4, cfgDec, numIter); 
        
        err1_dual = err1_dual + sum(decoded1(1:length(bits1)) ~= bits1);
        err2_dual = err2_dual + sum(decoded2(1:length(bits2)) ~= bits2);
        err3_dual = err3_dual + sum(decoded3(1:length(bits3)) ~= bits3);
        err4_dual = err4_dual + sum(decoded4(1:length(bits4)) ~= bits4); 
    end
    
    ber_dual_1(snrIdx) = err1_dual / (numTrials * infoLength);
    ber_dual_2(snrIdx) = err2_dual / (numTrials * infoLength);
    ber_dual_3(snrIdx) = err3_dual / (numTrials * infoLength);
    ber_dual_4(snrIdx) = err4_dual / (numTrials * infoLength); 
    ber_dual_total(snrIdx) = (err1_dual + err2_dual + err3_dual + err4_dual) / (numTrials * infoLength * streams); 
    
    fprintf('SNR=%2d: 交织LDPC S1=%.2e S2=%.2e S3=%.2e S4=%.2e | 独立LDPC S1=%.2e S2=%.2e S3=%.2e S4=%.2e\n', ...
        snr, ber_interleaved_1(snrIdx), ber_interleaved_2(snrIdx), ber_interleaved_3(snrIdx), ber_interleaved_4(snrIdx), ...
        ber_dual_1(snrIdx), ber_dual_2(snrIdx), ber_dual_3(snrIdx), ber_dual_4(snrIdx));
end

%画图
figure('Position', [100, 100, 1200, 600]);

subplot(1,2,1);
semilogy(snrRange, ber_interleaved_total, 'b-o', 'LineWidth', 2, 'DisplayName','Overall BER');
hold on; 
semilogy(snrRange, ber_interleaved_1, 'c--^', 'LineWidth', 1.5, 'DisplayName','Stream 1 (H=1)');
semilogy(snrRange, ber_interleaved_2, 'm--v', 'LineWidth', 1.5, 'DisplayName','Stream 2 (H=0.9)'); 
semilogy(snrRange, ber_interleaved_3, 'g--x', 'LineWidth', 1.5, 'DisplayName','Stream 3 (H=0.5)'); 
semilogy(snrRange, ber_interleaved_4, 'k--d', 'LineWidth', 1.5, 'DisplayName','Stream 4 (H=0.2)'); 
grid on; 
xlabel('Eb/N0 (dB)'); ylabel('Bit Error Rate (BER)'); 
title({'Architecture 1: Interleaved LDPC (4-Stream)'}); 
legend('show', 'Location', 'southwest'); 
ylim([1e-6 1]); 

subplot(1,2,2);
semilogy(snrRange, ber_dual_total, 'r-s', 'LineWidth', 2, 'DisplayName','Overall BER');
hold on; 
semilogy(snrRange, ber_dual_1, 'y-d', 'LineWidth', 1.5, 'DisplayName','Stream 1 (H=1)');
semilogy(snrRange, ber_dual_2, 'k-p', 'LineWidth', 1.5, 'DisplayName','Stream 2 (H=0.9)'); 
semilogy(snrRange, ber_dual_3, 'b-h', 'LineWidth', 1.5, 'DisplayName','Stream 3 (H=0.5)'); 
semilogy(snrRange, ber_dual_4, 'm-+', 'LineWidth', 1.5, 'DisplayName','Stream 4 (H=0.2)'); 
grid on;
xlabel('Eb/N0 (dB)'); ylabel('Bit Error Rate (BER)');
title({'Architecture 2: Original Dual LDPC (4-Stream)'});
legend('show', 'Location', 'southwest');
ylim([1e-6 1]);
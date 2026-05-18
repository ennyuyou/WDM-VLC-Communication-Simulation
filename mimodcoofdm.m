clc; clear; close all;

%% ============================================================
%  Interleaved-LDPC + 3×3 MIMO + DCO-OFDM (N=2048, 16QAM)
%  - 频域：使用正频 1..N/2-1 全部数据子载波（去 DC 与 Nyquist）
%  - 时域：IFFT 后加直流偏置 B，再做非负截断（DCO）
%  - 接收：先去偏置（减均值），再 FFT、等化、软解、LDPC 解码
%  - 只保留“交织 LDPC”路线
% ============================================================

%% ---------------- 可调参数 ----------------
M = 16;                          % QAM 阶数
bitsPerSymbol = log2(M);         % 每符号比特数 (16QAM=4)

N = 2048;                        % OFDM 点数（偶数）
data_carriers = 1:(N/2 - 1);     % DCO：用掉除 DC(0) 与 Nyquist(N/2) 外的全部正频
numSubcarriers = numel(data_carriers);   % = 1023

streams = 3;                     % 数据流 / 通道数
numTrials = 10000;               % 每个 SNR 仿真帧数（调试时可先小一点）
snrRange = 0:0.5:14;             % Eb/N0 (dB)

% 3x3 对角 MIMO 信道（互不串扰；仅增益不同）
H = diag([1, 0.8, 0.6]);

% DCO 偏置强度（经验 2.5~3.5；越大越少剪裁，但光功率开销更大）
kBias = 3.0;

% （可选）发射端 RMS 归一开关：对比 ACO/DCO 的电功率公平性
DO_RMS_NORMALIZE = false;

%% ---------------- LDPC 设置（802.11n 原型示例） ----------------
% 码长=648, 信息位=324, 码率=0.5
[cfgEnc, cfgDec] = getProtoMatrix(648, 324);
infoLength  = cfgEnc.NumInformationBits;       % 324
blockLength = cfgEnc.BlockLength;              % 648
numIter = 25;                                  % LDPC 最大迭代

fprintf('DCO-OFDM | 使用子载波数: %d | 数据流: %d\n', numSubcarriers, streams);

%% ---------------- 结果数组（仅交织 LDPC） ----------------
ber_interleaved_1     = zeros(1, length(snrRange));
ber_interleaved_2     = zeros(1, length(snrRange));
ber_interleaved_3     = zeros(1, length(snrRange));
ber_interleaved_total = zeros(1, length(snrRange));

%% ========= 在 parfor 外声明本地函数，以便并行调用 =========
% 放在脚本底部也可以；此处提前声明函数原型便于阅读
% function merged = interleave_bits(a,b,c)
% function [a,b,c] = deinterleave_bits(b_total, streams)

%% ---------------- 并行仿真（交织 LDPC） ----------------
parfor snrIdx = 1:length(snrRange)
    snr = snrRange(snrIdx);              % Eb/N0 (dB)
    EbN0_lin = 10^(snr/10);

    % Es/N0 = Eb/N0 * R * k
    codeRate = infoLength / blockLength;         % =0.5
    EsN0_lin = EbN0_lin * codeRate * bitsPerSymbol;

    % 噪声方差：复基带通常用 1/(2*EsN0)；DCO 我们在时域加“实噪声”
    % 这样定义后，FFT 到频域的符号噪声标度与之前流程兼容
    noiseVar = 1/(2*EsN0_lin);

    % 三路错误计数
    err1 = 0;  err2 = 0;  err3 = 0;

    for trial = 1:numTrials
        %% 1) 生成三路信息位（均分）
        % 324/3 = 108 比特每路
        data1 = randi([0,1], infoLength/streams, 1);
        data2 = randi([0,1], infoLength/streams, 1);
        data3 = randi([0,1], infoLength/streams, 1);

        %% 2) 交织 → LDPC 编码 → 解交织（恢复三路码字）
        merged   = interleave_bits(data1, data2, data3);     % [324 x 1]
        codeword = ldpcEncode(merged, cfgEnc);               % [648 x 1]
        [cw1, cw2, cw3] = deinterleave_bits(codeword, streams);

        %% 3) 按调制位数补齐（每路）
        pad = mod(-length(cw1), bitsPerSymbol);
        if pad ~= 0
            z = zeros(pad,1);
            cw1 = [cw1; z];  cw2 = [cw2; z];  cw3 = [cw3; z];
        end

        %% 4) 16QAM 调制（UnitAveragePower=true）
        s1 = qammod(cw1, M, 'InputType','bit','UnitAveragePower',true);
        s2 = qammod(cw2, M, 'InputType','bit','UnitAveragePower',true);
        s3 = qammod(cw3, M, 'InputType','bit','UnitAveragePower',true);

        %% 5) 限制每路符号不超过可用子载波（DCO: 1023）
        K = min([length(s1), length(s2), length(s3), numSubcarriers]);
        s1 = s1(1:K);  s2 = s2(1:K);  s3 = s3(1:K);
        data_use = data_carriers(1:K);
        idx = data_use + 1;                    % MATLAB 1 基索引

        %% 6) DCO-OFDM 频域映射 + IFFT + 偏置 + 非负截断
        % 频域厄米对称，确保时域实数
        X1 = zeros(N,1);  X2 = zeros(N,1);  X3 = zeros(N,1);
        X1(idx) = s1;     X1(N - data_use + 1) = conj(s1);
        X2(idx) = s2;     X2(N - data_use + 1) = conj(s2);
        X3(idx) = s3;     X3(N - data_use + 1) = conj(s3);

        x1 = real(ifft(X1) * sqrt(N));
        x2 = real(ifft(X2) * sqrt(N));
        x3 = real(ifft(X3) * sqrt(N));

        % 偏置：B = kBias * sigma（各路分别估计）
        B1 = kBias * std(x1);  B2 = kBias * std(x2);  B3 = kBias * std(x3);

        % 非负化（剪裁）
        x1_dco = max(x1 + B1, 0);
        x2_dco = max(x2 + B2, 0);
        x3_dco = max(x3 + B3, 0);

        % （可选）RMS 归一，保证与未偏置时的电功率相近（公平比较）
        if DO_RMS_NORMALIZE
            g1 = rms(x1) / max(rms(x1_dco), eps);
            g2 = rms(x2) / max(rms(x2_dco), eps);
            g3 = rms(x3) / max(rms(x3_dco), eps);
            x1_dco = g1 * x1_dco;  x2_dco = g2 * x2_dco;  x3_dco = g3 * x3_dco;
        end

        %% 7) MIMO 信道 + 实噪声
        tx = [x1_dco, x2_dco, x3_dco];            % N×3（非负实数）
        n  = sqrt(noiseVar) * randn(N, streams);   % 实噪声（更贴近 IM/DD）
        rx = tx * H + n;                           % N×3

        %% 8) 去偏置（均值法），FFT，取子载波并等化
        % 偏置通过信道变成 H(i,i)*B_i，均值法相当于估计并去除该分量
        rx1 = rx(:,1) - mean(rx(:,1));
        rx2 = rx(:,2) - mean(rx(:,2));
        rx3 = rx(:,3) - mean(rx(:,3));

        Y1 = fft(rx1) / sqrt(N);
        Y2 = fft(rx2) / sqrt(N);
        Y3 = fft(rx3) / sqrt(N);

        r1 = Y1(idx) / H(1,1);
        r2 = Y2(idx) / H(2,2);
        r3 = Y3(idx) / H(3,3);

        %% 9) 软解调（LLR）
        % 等化后噪声方差缩放
        effN1 = noiseVar / (abs(H(1,1))^2);
        effN2 = noiseVar / (abs(H(2,2))^2);
        effN3 = noiseVar / (abs(H(3,3))^2);

        llr1 = qamdemod(r1, M, 'OutputType','approxllr','UnitAveragePower',true,'NoiseVariance',effN1);
        llr2 = qamdemod(r2, M, 'OutputType','approxllr','UnitAveragePower',true,'NoiseVariance',effN2);
        llr3 = qamdemod(r3, M, 'OutputType','approxllr','UnitAveragePower',true,'NoiseVariance',effN3);

        % 截回调制前长度（去掉 padding 部分）
        llr1 = llr1(1:length(cw1));
        llr2 = llr2(1:length(cw2));
        llr3 = llr3(1:length(cw3));

        %% 10) LLR 交织 → LDPC 解码 → 解交织
        llr_total     = interleave_bits(llr1, llr2, llr3);
        decoded_total = ldpcDecode(llr_total, cfgDec, numIter);
        [dec1, dec2, dec3] = deinterleave_bits(decoded_total, streams);

        %% 11) 与原始信息位比较（只比 info 位长度）
        err1 = err1 + sum(dec1(1:length(data1)) ~= data1);
        err2 = err2 + sum(dec2(1:length(data2)) ~= data2);
        err3 = err3 + sum(dec3(1:length(data3)) ~= data3);
    end

    %% 12) BER 统计（交织 LDPC）
    ber_interleaved_1(snrIdx)     = err1 / (numTrials * infoLength/streams);
    ber_interleaved_2(snrIdx)     = err2 / (numTrials * infoLength/streams);
    ber_interleaved_3(snrIdx)     = err3 / (numTrials * infoLength/streams);
    ber_interleaved_total(snrIdx) = (err1 + err2 + err3) / (numTrials * infoLength);

    fprintf('SNR=%4.1f dB | Interleaved-LDPC (DCO): S1=%.2e  S2=%.2e  S3=%.2e  | Total=%.2e\n', ...
        snr, ber_interleaved_1(snrIdx), ber_interleaved_2(snrIdx), ber_interleaved_3(snrIdx), ber_interleaved_total(snrIdx));
end

%% ---------------- 画图（仅交织 LDPC） ----------------
figure('Position',[100,100,900,500]);
semilogy(snrRange, ber_interleaved_total, 'b-o', 'LineWidth', 2, 'DisplayName','Overall BER');
hold on; grid on;
semilogy(snrRange, ber_interleaved_1, 'c--^', 'LineWidth', 1.5, 'DisplayName','Stream 1 (H=1)');
semilogy(snrRange, ber_interleaved_2, 'm--v', 'LineWidth', 1.5, 'DisplayName','Stream 2 (H=0.8)');
semilogy(snrRange, ber_interleaved_3, 'g--x', 'LineWidth', 1.5, 'DisplayName','Stream 3 (H=0.6)');
xlabel('Eb/N0 (dB)'); ylabel('Bit Error Rate (BER)');
title('Interleaved LDPC + DCO-OFDM (3 Streams)');
legend('show','Location','southwest');
ylim([1e-6 1]); xlim([snrRange(1) snrRange(end)]);

%% ================== 本地函数 ==================
function merged = interleave_bits(a, b, c)
% a,b,c 交织成 a1,b1,c1,a2,b2,c2,...
    merged = reshape([a(:)'; b(:)'; c(:)'], [], 1);
end

function [a, b, c] = deinterleave_bits(b_total, streams)
% 从交织序列复原三路（streams=3）
    a = b_total(1:streams:end);
    b = b_total(2:streams:end);
    c = b_total(3:streams:end);
end

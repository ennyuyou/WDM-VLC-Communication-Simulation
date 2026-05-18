% DCO-OFDM Receiver WDM interleaving
% Author: Akatsuki Sky
clc; 
clear all;
close all;

%% ---------------- 路径与参数加载 ----------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

allPath = fullfile(scriptDir, 'Knowledge.mat');

load(allPath, ...
    'numChannels','M','bitsPerSymbol','N','Ncp','numsym_data','numsym_pilot', ...
    'pilotLen','dataLen', ...
    'finalSerialSignal','finalSerialSignal_ch2','finalSerialSignal_ch3', ...
    'upsampleRate','s_vec_tx_pilot','posLower', ...
    'K','blockLength','numIter', ...
    'txBits1','txBits2','txBits3','numBitsPerTX');


numDataSubcarriers = numel(posLower);


% orgrx1 = readmatrix('tx1_simple.csv');   % RX1
% orgrx2 = readmatrix('tx2_simple.csv');   % RX2
% orgrx3 = readmatrix('tx3_simple.csv');   % RX3


% RXab means from transmitter a to receiver b

RX11 = readmatrix('RigolDS6.csv'); 
RX12 = readmatrix('RigolDS7.csv'); 
RX13 = readmatrix('RigolDS15.csv');

RX21 = readmatrix('RigolDS8.csv');
RX22 = readmatrix('RigolDS9.csv'); 
RX23 = readmatrix('RigolDS16.csv');

RX31 = readmatrix('RigolDS10.csv'); 
RX32 = readmatrix('RigolDS11.csv'); 
RX33 = readmatrix('RigolDS17.csv'); 



% RXCTa means crosstalk complex on receiver a
upsampleRate = 5;

[RX11idx, RX11sync] = syncWithReference(RX11, finalSerialSignal, upsampleRate);
[RX12idx, RX12sync] = syncWithReference(RX12, finalSerialSignal, upsampleRate);
[RX13idx, RX13sync] = syncWithReference(RX13, finalSerialSignal, upsampleRate);

[RX21idx, RX21sync] = syncWithReference(RX21, finalSerialSignal_ch2, upsampleRate);
[RX22idx, RX22sync] = syncWithReference(RX22, finalSerialSignal_ch2, upsampleRate);
[RX23idx, RX23sync] = syncWithReference(RX23, finalSerialSignal_ch2, upsampleRate);

[RX31idx, RX31sync] = syncWithReference(RX31, finalSerialSignal_ch3, upsampleRate);
[RX32idx, RX32sync] = syncWithReference(RX32, finalSerialSignal_ch3, upsampleRate);
[RX33idx, RX33sync] = syncWithReference(RX33, finalSerialSignal_ch3, upsampleRate);

PD1=RX11sync+RX21sync+RX31sync;
PD2=RX12sync+RX22sync+RX32sync;
PD3=RX13sync+RX23sync+RX33sync;

% PD1=RX11sync;
% PD2=RX22sync;
% PD3=RX33sync;


rx1_sync = PD1;
rx2_sync = PD2;
rx3_sync = PD3;

% %% ================== 添加人为高斯白噪声 ==================
% % 设为 10, 15, 20 等数值添加噪声；设为 inf (无穷大) 则不加噪声。
% add_SNR_dB = 10; 
% 
% if ~isinf(add_SNR_dB)
%     sigPower1 = mean(rx1_sync.^2);
%     sigPower2 = mean(rx2_sync.^2);
%     sigPower3 = mean(rx3_sync.^2);
%     noisePower1 = sigPower1 / (10^(add_SNR_dB/10));
%     noisePower2 = sigPower2 / (10^(add_SNR_dB/10));
%     noisePower3 = sigPower3 / (10^(add_SNR_dB/10));
%     rx1_sync = rx1_sync + sqrt(noisePower1) * randn(size(rx1_sync));
%     rx2_sync = rx2_sync + sqrt(noisePower2) * randn(size(rx2_sync));
%     rx3_sync = rx3_sync + sqrt(noisePower3) * randn(size(rx3_sync));
% 
%     fprintf('已添加人为噪声: Target SNR = %.1f dB\n', add_SNR_dB);
% end
% %% ================== 不需要加噪声的时候必须注释掉 ==================

pilotRegionLen = 3 * pilotLen;

% ===== RX1 =====
pilot_rx1_slot1_time = rx1_sync(1                : pilotLen);
pilot_rx1_slot2_time = rx1_sync(pilotLen+1       : 2*pilotLen);
pilot_rx1_slot3_time = rx1_sync(2*pilotLen+1     : 3*pilotLen);
data_rx1_time        = rx1_sync(pilotRegionLen+1 : pilotRegionLen+dataLen);
            
% ===== RX2 =====
pilot_rx2_slot1_time = rx2_sync(1                : pilotLen);
pilot_rx2_slot2_time = rx2_sync(pilotLen+1       : 2*pilotLen);
pilot_rx2_slot3_time = rx2_sync(2*pilotLen+1     : 3*pilotLen);
data_rx2_time        = rx2_sync(pilotRegionLen+1 : pilotRegionLen+dataLen);

% ===== RX3 =====
pilot_rx3_slot1_time = rx3_sync(1                : pilotLen);
pilot_rx3_slot2_time = rx3_sync(pilotLen+1       : 2*pilotLen);
pilot_rx3_slot3_time = rx3_sync(2*pilotLen+1     : 3*pilotLen);
data_rx3_time        = rx3_sync(pilotRegionLen+1 : pilotRegionLen+dataLen);

%% ---------------- 发送端导频：时域 -> 频域 ----------------
% s_vec_tx_pilot 是发射端导频（所有符号连成一串，带 CP）
pilot_tx_time_grid = reshape(s_vec_tx_pilot, N + Ncp, numsym_pilot);  
pilot_tx_time_noCP = pilot_tx_time_grid(Ncp+1:end, :);              
pilot_TX_freq_full = fft(pilot_tx_time_noCP, N, 1) / sqrt(N);       
pilot_TX_freq_sub  = pilot_TX_freq_full(posLower, :);             


% ===== RX1：三个导频的频域（子载波维度） =====
pilot_rx1_slot1_grid = reshape(pilot_rx1_slot1_time, N + Ncp, numsym_pilot);
pilot_rx1_slot1_noCP = pilot_rx1_slot1_grid(Ncp+1:end, :);
pilot_RX1_slot1_full = fft(pilot_rx1_slot1_noCP, N, 1) / sqrt(N);
pilot_RX1_1  = pilot_RX1_slot1_full(posLower, :);   

pilot_rx1_slot2_grid = reshape(pilot_rx1_slot2_time, N + Ncp, numsym_pilot);
pilot_rx1_slot2_noCP = pilot_rx1_slot2_grid(Ncp+1:end, :);
pilot_RX1_slot2_full = fft(pilot_rx1_slot2_noCP, N, 1) / sqrt(N);
pilot_RX1_2  = pilot_RX1_slot2_full(posLower, :);

pilot_rx1_slot3_grid = reshape(pilot_rx1_slot3_time, N + Ncp, numsym_pilot);
pilot_rx1_slot3_noCP = pilot_rx1_slot3_grid(Ncp+1:end, :);
pilot_RX1_slot3_full = fft(pilot_rx1_slot3_noCP, N, 1) / sqrt(N);
pilot_RX1_3  = pilot_RX1_slot3_full(posLower, :);

% ===== RX2 =====
pilot_rx2_slot1_grid = reshape(pilot_rx2_slot1_time, N + Ncp, numsym_pilot);
pilot_rx2_slot1_noCP = pilot_rx2_slot1_grid(Ncp+1:end, :);
pilot_RX2_slot1_full = fft(pilot_rx2_slot1_noCP, N, 1) / sqrt(N);
pilot_RX2_1  = pilot_RX2_slot1_full(posLower, :);

pilot_rx2_slot2_grid = reshape(pilot_rx2_slot2_time, N + Ncp, numsym_pilot);
pilot_rx2_slot2_noCP = pilot_rx2_slot2_grid(Ncp+1:end, :);
pilot_RX2_slot2_full = fft(pilot_rx2_slot2_noCP, N, 1) / sqrt(N);
pilot_RX2_2  = pilot_RX2_slot2_full(posLower, :);

pilot_rx2_slot3_grid = reshape(pilot_rx2_slot3_time, N + Ncp, numsym_pilot);
pilot_rx2_slot3_noCP = pilot_rx2_slot3_grid(Ncp+1:end, :);
pilot_RX2_slot3_full = fft(pilot_rx2_slot3_noCP, N, 1) / sqrt(N);
pilot_RX2_3  = pilot_RX2_slot3_full(posLower, :);

% ===== RX3 =====
pilot_rx3_slot1_grid = reshape(pilot_rx3_slot1_time, N + Ncp, numsym_pilot);
pilot_rx3_slot1_noCP = pilot_rx3_slot1_grid(Ncp+1:end, :);
pilot_RX3_slot1_full = fft(pilot_rx3_slot1_noCP, N, 1) / sqrt(N);
pilot_RX3_1  = pilot_RX3_slot1_full(posLower, :);

pilot_rx3_slot2_grid = reshape(pilot_rx3_slot2_time, N + Ncp, numsym_pilot);
pilot_rx3_slot2_noCP = pilot_rx3_slot2_grid(Ncp+1:end, :);
pilot_RX3_slot2_full = fft(pilot_rx3_slot2_noCP, N, 1) / sqrt(N);
pilot_RX3_2  = pilot_RX3_slot2_full(posLower, :);

pilot_rx3_slot3_grid = reshape(pilot_rx3_slot3_time, N + Ncp, numsym_pilot);
pilot_rx3_slot3_noCP = pilot_rx3_slot3_grid(Ncp+1:end, :);
pilot_RX3_slot3_full = fft(pilot_rx3_slot3_noCP, N, 1) / sqrt(N);
pilot_RX3_3  = pilot_RX3_slot3_full(posLower, :);


%% ---------------- 时域 -> 频域（子载波） ----------------

% RX1 数据
data_rx1_grid = reshape(data_rx1_time, N + Ncp, numsym_data);
data_rx1_noCP = data_rx1_grid(Ncp+1:end, :);
data_rx1_full = fft(data_rx1_noCP, N, 1) / sqrt(N);
data_rx1  = data_rx1_full(posLower, :);  

% RX2 数据
data_rx2_grid = reshape(data_rx2_time, N + Ncp, numsym_data);
data_rx2_noCP = data_rx2_grid(Ncp+1:end, :);
data_rx2_full = fft(data_rx2_noCP, N, 1) / sqrt(N);
data_rx2  = data_rx2_full(posLower, :);

% RX3 数据
data_rx3_grid = reshape(data_rx3_time, N + Ncp, numsym_data);
data_rx3_noCP = data_rx3_grid(Ncp+1:end, :);
data_rx3_full = fft(data_rx3_noCP, N, 1) / sqrt(N);
data_rx3  = data_rx3_full(posLower, :);




figure;
title(' TX1 星座 均衡前');
dscatter(reshape(real(data_rx1),[],1), reshape(imag(data_rx1),[],1));

figure;
title(' TX2 星座 均衡前');
dscatter(reshape(real(data_rx2),[],1), reshape(imag(data_rx2),[],1));

figure;
title(' TX3 星座 均衡前');
dscatter(reshape(real(data_rx3),[],1), reshape(imag(data_rx3),[],1));



eps_val = 1e-12;

% 先做"相除"：R_ij(k,n) = RX_ij(k,n) / TX(k,n)
R11 = pilot_RX1_1 ./ (pilot_TX_freq_sub + eps_val);
R12 = pilot_RX1_2 ./ (pilot_TX_freq_sub + eps_val);
R13 = pilot_RX1_3 ./ (pilot_TX_freq_sub + eps_val);

R21 = pilot_RX2_1 ./ (pilot_TX_freq_sub + eps_val);
R22 = pilot_RX2_2 ./ (pilot_TX_freq_sub + eps_val);
R23 = pilot_RX2_3 ./ (pilot_TX_freq_sub + eps_val);

R31 = pilot_RX3_1 ./ (pilot_TX_freq_sub + eps_val);
R32 = pilot_RX3_2 ./ (pilot_TX_freq_sub + eps_val);
R33 = pilot_RX3_3 ./ (pilot_TX_freq_sub + eps_val);

% 对导频符号维度（第 2 维）取平均：
% H_ij(k) = mean_n( R_ij(k,n) )
H11 = mean(R11, 2);    % numDataSubcarriers × 1
H12 = mean(R12, 2);
H13 = mean(R13, 2);

H21 = mean(R21, 2);
H22 = mean(R22, 2);
H23 = mean(R23, 2);

H31 = mean(R31, 2);
H32 = mean(R32, 2);
H33 = mean(R33, 2);

H_all = [H11 H12 H13 H21 H22 H23 H31 H32 H33];

numSC = numDataSubcarriers;

% size(H) = 3 × 3 × 127
% premute: 重新排列矩阵
H = permute(reshape(H_all, numSC, 3, 3), [2 3 1]);


avgH = squeeze(mean(H, 3)).';
fprintf('H 信道估计平均值 (行=Rx, 列=Tx):\n');
disp(abs(avgH));
fprintf('svd(H) :\n');
svd(avgH)



numSC  = numDataSubcarriers;
numsym = numsym_data;

% 预分配：numSC × numsym × 3
% size:  127 × 200 × 3
Xhat = zeros(numSC, numsym, 3);

for k = 1:numSC
    % size(Yk) = 200,3
    % 代表位于当前子载波的三个信道的数据
    Yk = [data_rx1(k, :); 
          data_rx2(k, :); 
          data_rx3(k, :)].';

    Xk = Yk / H(:, :, k);
    Xhat(k, :, :) = Xk;
end

eq_data1 = squeeze(Xhat(:, :, 1)); 
eq_data2 = squeeze(Xhat(:, :, 2));
eq_data3 = squeeze(Xhat(:, :, 3));


figure;
title(' TX1 星座 均衡后');
dscatter(reshape(real(eq_data1),[],1), reshape(imag(eq_data1),[],1));
axis([-4 4 -4 4]);
figure;
title(' TX2 星座 均衡后');
dscatter(reshape(real(eq_data2),[],1), reshape(imag(eq_data2),[],1));
axis([-4 4 -4 4]);
figure;
title(' TX3 星座 均衡后');
dscatter(reshape(real(eq_data3),[],1), reshape(imag(eq_data3),[],1));
axis([-4 4 -4 4]);


numSC_p = numDataSubcarriers;
Np      = numsym_pilot;

rxeq1 = zeros(numSC_p, Np);  
rxeq2 = zeros(numSC_p, Np);  
rxeq3 = zeros(numSC_p, Np);   

for k = 1:numSC_p
    Hk = H(:, :, k);      % 3×3

    Y1 = [ pilot_RX1_1(k, :);
           pilot_RX2_1(k, :);
           pilot_RX3_1(k, :) ].';      % Np×3
    X1_hat = Y1 / Hk;                  % Np×3，等价于 Y1 * inv(Hk)
    rxeq1(k, :) = X1_hat(:, 1).';      % 第 1 列对应 TX1

    Y2 = [ pilot_RX1_2(k, :);
           pilot_RX2_2(k, :);
           pilot_RX3_2(k, :) ].';
    X2_hat = Y2 / Hk;                  % Np×3
    rxeq2(k, :) = X2_hat(:, 2).';      % 第 2 列对应 TX2

    Y3 = [ pilot_RX1_3(k, :);
           pilot_RX2_3(k, :);
           pilot_RX3_3(k, :) ].';
    X3_hat = Y3 / Hk;                  % Np×3
    rxeq3(k, :) = X3_hat(:, 3).';      % 第 3 列对应 TX3
end

% 做差
diff1 = rxeq1 - pilot_TX_freq_sub;  
diff2 = rxeq2 - pilot_TX_freq_sub; 
diff3 = rxeq3 - pilot_TX_freq_sub;

% 算成一个 1×3 向量
noiseVar_raw = [
    mean(abs(diff1(:)).^2, 'omitnan'), ...
    mean(abs(diff2(:)).^2, 'omitnan'), ...
    mean(abs(diff3(:)).^2, 'omitnan')
];

pilotPower = mean(abs(pilot_TX_freq_sub(:)).^2, 'omitnan'); 
noiseFloor = 1e-2 * pilotPower; 
noiseVar = max(noiseVar_raw, noiseFloor);
noiseVar(~isfinite(noiseVar)) = noiseFloor;


noiseVar1 = noiseVar(1);
noiseVar2 = noiseVar(2);
noiseVar3 = noiseVar(3);



SNR_est_1 = 10 * log10(pilotPower / noiseVar1);
SNR_est_2 = 10 * log10(pilotPower / noiseVar2);
SNR_est_3 = 10 * log10(pilotPower / noiseVar3);

fprintf('------------------------------------------------\n');
fprintf('Estimated SNR (based on Pilots):\n');
fprintf('Channel 1: %.2f dB\n', SNR_est_1);
fprintf('Channel 2: %.2f dB\n', SNR_est_2);
fprintf('Channel 3: %.2f dB\n', SNR_est_3);
fprintf('------------------------------------------------\n');



llr1 = qamdemod(eq_data1(:), M, 'gray', ...
    'OutputType','llr', ...
    'UnitAveragePower', false, ...
    'NoiseVariance', noiseVar1);

llr2 = qamdemod(eq_data2(:), M, 'gray', ...
    'OutputType','llr', ...
    'UnitAveragePower', false, ...
    'NoiseVariance', noiseVar2);

llr3 = qamdemod(eq_data3(:), M, 'gray', ...
    'OutputType','llr', ...
    'UnitAveragePower', false, ...
    'NoiseVariance', noiseVar3);


%% ================== C. 联合 LLR 交织与 LDPC 解码 ==================
[cfgEnc_rx, cfgDec_rx] = getProtoMatrix(648, 540);

% 1. 将每路的 LLR 重塑为 [bitsPerSymbol, numSymbols]
L1 = reshape(llr1, bitsPerSymbol, []);
L2 = reshape(llr2, bitsPerSymbol, []);
L3 = reshape(llr3, bitsPerSymbol, []);

% 2. 垂直堆叠：形成 [Ch1_Sym1; Ch2_Sym1; Ch3_Sym1; ...] 的结构
% 结果矩阵的行数为 3 * bitsPerSymbol
L_stack = [L1; L2; L3];

% 3. 按列展开，恢复成发射端编码后的比特流顺序
llr_combined = L_stack(:);

% 插入全局解交织 (De-interleaving)
rng(42); % 1. 强制重置种子 (必须与 Tx 一致！)

% 计算符号总数
len_syms = length(llr_combined) / bitsPerSymbol;
perm_idx = randperm(len_syms);

% 2. 把 LLR 变回 [bitsPerSymbol x NumSyms] 矩阵以进行符号级操作
llr_per_symbol_matrix = reshape(llr_combined, bitsPerSymbol, []);

% 3. 计算逆映射索引
inv_perm_idx = zeros(1, len_syms);
inv_perm_idx(perm_idx) = 1:len_syms;

% 4. 执行解交织 (列交换)
llr_restored_matrix = llr_per_symbol_matrix(:, inv_perm_idx);

% 5. 变回比特流
llr_final_stream = llr_restored_matrix(:);

% 6. 解码 (【注意】输入必须是 llr_final_stream)
[decBits_combined, ~] = decodeLdpcStreamSimple(llr_final_stream, K, blockLength, cfgDec_rx, numIter);

totalInfoBits = length(txBits1) + length(txBits2) + length(txBits3);
% 截取有效数据
if length(decBits_combined) < totalInfoBits
    warning('解码输出长度 (%d) 小于预期信息位 (%d)，可能误码率极高导致丢包', length(decBits_combined), totalInfoBits);
    % 补零以防报错，方便看 BER
    decBits_combined = [decBits_combined; zeros(totalInfoBits - length(decBits_combined), 1)];
end
decBits_combined = decBits_combined(1:totalInfoBits);

% 拆分回三路
bits_reshaped = reshape(decBits_combined, 3, []); % 3 x N
decBits1 = bits_reshaped(1, :).';
decBits2 = bits_reshaped(2, :).';
decBits3 = bits_reshaped(3, :).';

%% ================== D. 计算 BER ==================
len1 = length(decBits1);
len2 = length(decBits2);
len3 = length(decBits3);

ref1 = txBits1(1:len1);
ref2 = txBits2(1:len2);
ref3 = txBits3(1:len3);

err1 = nnz(ref1 ~= decBits1);
err2 = nnz(ref2 ~= decBits2);
err3 = nnz(ref3 ~= decBits3);

BER1 = err1 / max(len1, 1);
BER2 = err2 / max(len2, 1);
BER3 = err3 / max(len3, 1);

BER_total = (err1 + err2 + err3) / max(len1 + len2 + len3, 1);

fprintf('BER_TX1 = %.3e  (errors = %d / %d bits)\n', BER1, err1, len1);
fprintf('BER_TX2 = %.3e  (errors = %d / %d bits)\n', BER2, err2, len2);
fprintf('BER_TX3 = %.3e  (errors = %d / %d bits)\n', BER3, err3, len3);
fprintf('BER_Total = %.3e\n', BER_total);




function [syncIndex, synced] = syncWithReference(receivedSignal, referenceSignal, upsampleRate)
    original = repelem(referenceSignal, upsampleRate);
    [c, lags] = xcorr(receivedSignal, original);
    [~, idx]  = max(c);


    syncIndex = lags(idx) + 1;
    len = length(referenceSignal) * upsampleRate;
    syncIndex = max(syncIndex, 1);
    maxStart  = numel(receivedSignal) - len + 1;
    syncIndex = min(syncIndex, maxStart);

    if syncIndex < 1 || syncIndex + len - 1 > numel(receivedSignal)
        error('同步失败');
    end
    syncedSegment = receivedSignal(syncIndex : syncIndex + len - 1);
    synced = syncedSegment(3:upsampleRate:end);
end

function [infoBitsVec, numBlocksDecoded] = decodeLdpcStreamSimple(llrStream, K, blockLength, cfgDec, numIter)

    LLR_CLIP = 20;
    llrStream = max(min(llrStream, LLR_CLIP), -LLR_CLIP);

    numBlocksDecoded = floor(numel(llrStream) / blockLength);
    if numBlocksDecoded <= 0
        infoBitsVec = zeros(0,1);
        return;
    end

    llrStream = llrStream(1:numBlocksDecoded * blockLength);
    llrMatrix = reshape(llrStream, blockLength, numBlocksDecoded);

    infoBitsHat = zeros(K, numBlocksDecoded);
    for blk = 1:numBlocksDecoded
        infoBitsHat(:, blk) = ldpcDecode(llrMatrix(:, blk), cfgDec, numIter);
    end

    infoBitsVec = infoBitsHat(:);
end

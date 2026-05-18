% WDM_BER_Sweep.m
% 用你的 WDM_TX/WDM_RX 流程，扫 SNR 并画 BER 曲线
clc; clear; close all;

%% 0) 先跑一次 TX，生成 Knowledge.mat（包含 bits / reference / pilot 等）
run('WDM_TX.m');

%% 1) 加载 TX 保存的参数
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
allPath = fullfile(scriptDir, 'Knowledge.mat');

load(allPath, ...
    'M','bitsPerSymbol','N','Ncp','numsym_data','numsym_pilot', ...
    'pilotLen','dataLen', ...
    'finalSerialSignal','finalSerialSignal_ch2','finalSerialSignal_ch3', ...
    's_vec_tx_pilot','posLower', ...
    'K','blockLength','numIter', ...
    'txBits1','txBits2','txBits3');

numDataSubcarriers = numel(posLower);
pilotRegionLen = 3 * pilotLen;

%% 2) 读取 Rigol 数据（你的原 RX 文件读法）
RX11 = readmatrix('RigolDS6.csv'); 
RX12 = readmatrix('RigolDS7.csv'); 
RX13 = readmatrix('RigolDS15.csv');

RX21 = readmatrix('RigolDS8.csv');
RX22 = readmatrix('RigolDS9.csv'); 
RX23 = readmatrix('RigolDS16.csv');

RX31 = readmatrix('RigolDS10.csv'); 
RX32 = readmatrix('RigolDS11.csv'); 
RX33 = readmatrix('RigolDS17.csv'); 

%% 3) 同步（只做一次！）得到三路 PD1/2/3（与你 RX 一致）
upsampleRate = 5;

[~, RX11sync] = syncWithReference(RX11, finalSerialSignal, upsampleRate);
[~, RX12sync] = syncWithReference(RX12, finalSerialSignal, upsampleRate);
[~, RX13sync] = syncWithReference(RX13, finalSerialSignal, upsampleRate);

[~, RX21sync] = syncWithReference(RX21, finalSerialSignal_ch2, upsampleRate);
[~, RX22sync] = syncWithReference(RX22, finalSerialSignal_ch2, upsampleRate);
[~, RX23sync] = syncWithReference(RX23, finalSerialSignal_ch2, upsampleRate);

[~, RX31sync] = syncWithReference(RX31, finalSerialSignal_ch3, upsampleRate);
[~, RX32sync] = syncWithReference(RX32, finalSerialSignal_ch3, upsampleRate);
[~, RX33sync] = syncWithReference(RX33, finalSerialSignal_ch3, upsampleRate);

PD1 = RX11sync + RX21sync + RX31sync;
PD2 = RX12sync + RX22sync + RX32sync;
PD3 = RX13sync + RX23sync + RX33sync;

rx1_base = PD1;
rx2_base = PD2;
rx3_base = PD3;

%% 4) 扫 SNR
SNRdB_vec = 0:1:20;     % 你按需要改
Nframes   = 20;         % 每个 SNR 重复次数（想快就调小）

BER_tx1 = zeros(size(SNRdB_vec));
BER_tx2 = zeros(size(SNRdB_vec));
BER_tx3 = zeros(size(SNRdB_vec));
BER_all = zeros(size(SNRdB_vec));

for ii = 1:numel(SNRdB_vec)
    snrdb = SNRdB_vec(ii);

    err = [0 0 0];
    tot = [0 0 0];

    for ff = 1:Nframes
        % 4.1 给同步后的波形加人为噪声（基于当前波形功率定义的 SNR）
        rx1 = add_awgn_real(rx1_base, snrdb);
        rx2 = add_awgn_real(rx2_base, snrdb);
        rx3 = add_awgn_real(rx3_base, snrdb);

        % 4.2 调用"封装后的 RX 解码函数"
        [e, t] = wdm_rx_decode_once( ...
            rx1, rx2, rx3, ...
            M, bitsPerSymbol, N, Ncp, ...
            numsym_data, numsym_pilot, ...
            pilotLen, dataLen, pilotRegionLen, ...
            posLower, s_vec_tx_pilot, ...
            K, blockLength, numIter, ...
            txBits1, txBits2, txBits3);

        err = err + e;
        tot = tot + t;
    end

    BER_tx1(ii) = err(1)/max(tot(1),1);
    BER_tx2(ii) = err(2)/max(tot(2),1);
    BER_tx3(ii) = err(3)/max(tot(3),1);
    BER_all(ii) = sum(err)/max(sum(tot),1);

    fprintf('SNR=%.1f dB: BER=[%.3e %.3e %.3e], overall=%.3e\n', ...
        snrdb, BER_tx1(ii), BER_tx2(ii), BER_tx3(ii), BER_all(ii));
end

%% 5) 画 BER 曲线
figure('Color','w'); grid on; hold on;
semilogy(SNRdB_vec, BER_tx1, '-o', 'LineWidth', 1);
semilogy(SNRdB_vec, BER_tx2, '-s', 'LineWidth', 1);
semilogy(SNRdB_vec, BER_tx3, '-^', 'LineWidth', 1);
semilogy(SNRdB_vec, BER_all, '-d', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BER');
title('WDM DCO-OFDM BER vs SNR (based on your TX/RX chain)');
legend('TX1','TX2','TX3','overall','Location','northeast');
ylim([1e-5 1]);

%% ====== 工具函数：和你 RX 里一致 ======
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

function y = add_awgn_real(x, SNRdB)
    sigPower = mean(x.^2);
    noisePower = sigPower / (10^(SNRdB/10));
    y = x + sqrt(noisePower) * randn(size(x));
end

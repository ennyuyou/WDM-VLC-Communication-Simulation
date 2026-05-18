% WDM_BER_Normal.m
% 目标：正常 BER（semilogy），不做 BER 地板替换，不用 ylim 限制
clc; clear; close all;

%% 0) 跑 TX（生成 Knowledge.mat）
run('WDM_TX.m');

%% 1) 载入 TX 参数
load('Knowledge.mat', ...
    'M','bitsPerSymbol','N','Ncp','numsym_data','numsym_pilot', ...
    'pilotLen','dataLen','s_vec_tx_pilot','posLower', ...
    'K','blockLength','numIter', ...
    'txBits1','txBits2','txBits3', ...
    'finalSerialSignal','finalSerialSignal_ch2','finalSerialSignal_ch3');

pilotRegionLen = 3 * pilotLen;

%% 2) 读 Rigol
RX11 = readmatrix('RigolDS6.csv'); 
RX12 = readmatrix('RigolDS7.csv'); 
RX13 = readmatrix('RigolDS15.csv');

RX21 = readmatrix('RigolDS8.csv');
RX22 = readmatrix('RigolDS9.csv'); 
RX23 = readmatrix('RigolDS16.csv');

RX31 = readmatrix('RigolDS10.csv'); 
RX32 = readmatrix('RigolDS11.csv'); 
RX33 = readmatrix('RigolDS17.csv'); 

%% 3) 同步（只做一次）
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

rx1_base = PD1; rx2_base = PD2; rx3_base = PD3;

%% 4) 扫 SNR
SNRdB_vec = 0:0.5:15;   % 建议扫宽一点才有 waterfall
Nframes   = 100;         % 越大越平滑

BER_tx1 = zeros(size(SNRdB_vec));
BER_tx2 = zeros(size(SNRdB_vec));
BER_tx3 = zeros(size(SNRdB_vec));
BER_all = zeros(size(SNRdB_vec));

for ii = 1:numel(SNRdB_vec)
    snrdb = SNRdB_vec(ii);

    err = [0 0 0];
    tot = [0 0 0];

    for ff = 1:Nframes
        rx1 = add_awgn_real(rx1_base, snrdb);
        rx2 = add_awgn_real(rx2_base, snrdb);
        rx3 = add_awgn_real(rx3_base, snrdb);

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

    fprintf('SNR=%.2f dB: BER=[%.3e %.3e %.3e], overall=%.3e\n', ...
        snrdb, BER_tx1(ii), BER_tx2(ii), BER_tx3(ii), BER_all(ii));
end

%% 5) 画"正常"对数 BER（不改 BER 数值）
% semilogy 不能画 0，所以我们只对"显示"做处理：0 → NaN（不改变原 BER 数组）
p1 = BER_tx1; p2 = BER_tx2; p3 = BER_tx3; pall = BER_all;
p1(p1==0)=NaN; p2(p2==0)=NaN; p3(p3==0)=NaN; pall(pall==0)=NaN;

figure('Color','w','Position',[100 100 850 550]); hold on; grid on;
semilogy(SNRdB_vec, p1, '-o', 'LineWidth', 1, 'MarkerSize', 5);
semilogy(SNRdB_vec, p2, '-s', 'LineWidth', 1, 'MarkerSize', 5);
semilogy(SNRdB_vec, p3, '-^', 'LineWidth', 1, 'MarkerSize', 5);
semilogy(SNRdB_vec, pall,'-d', 'LineWidth', 2, 'MarkerSize', 5);

title('3-Channel DCO-OFDM', 'FontSize', 14, 'FontWeight','normal');
xlabel('$E_s/N_0$ (dB)', 'Interpreter','latex', 'FontSize', 13);
ylabel('BER', 'FontSize', 12);
ax = gca;
ax.YScale = 'log';
ax.YTick = 10.^(-5:0);   % 例如显示 10^-5 到 10^0
ax.YTickLabel = arrayfun(@(x) sprintf('10^{%d}', x), -5:0, 'UniformOutput', false);

set(gca,'YMinorGrid','on','XMinorGrid','on','FontSize',11,'LineWidth',0.8);
legend('Channel1','Channel2','Channel3','overallBER', 'Location','northeast');
hold off;

%% ========= 工具函数 =========
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
        error('同步失败：syncIndex 越界');
    end

    syncedSegment = receivedSignal(syncIndex : syncIndex + len - 1);
    synced = syncedSegment(3:upsampleRate:end);
end

function y = add_awgn_real(x, SNRdB)
    sigPower   = mean(x.^2);
    noisePower = sigPower / (10^(SNRdB/10));
    y = x + sqrt(noisePower) * randn(size(x));
end

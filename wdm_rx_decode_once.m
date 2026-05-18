function [err3, tot3] = wdm_rx_decode_once( ...
    rx1_sync, rx2_sync, rx3_sync, ...
    M, bitsPerSymbol, N, Ncp, ...
    numsym_data, numsym_pilot, ...
    pilotLen, dataLen, pilotRegionLen, ...
    posLower, s_vec_tx_pilot, ...
    K, blockLength, numIter, ...
    txBits1, txBits2, txBits3)

% === 1) 切分导频/数据（与你 RX 一致） ===
pilot_rx1_slot1_time = rx1_sync(1                : pilotLen);
pilot_rx1_slot2_time = rx1_sync(pilotLen+1       : 2*pilotLen);
pilot_rx1_slot3_time = rx1_sync(2*pilotLen+1     : 3*pilotLen);
data_rx1_time        = rx1_sync(pilotRegionLen+1 : pilotRegionLen+dataLen);

pilot_rx2_slot1_time = rx2_sync(1                : pilotLen);
pilot_rx2_slot2_time = rx2_sync(pilotLen+1       : 2*pilotLen);
pilot_rx2_slot3_time = rx2_sync(2*pilotLen+1     : 3*pilotLen);
data_rx2_time        = rx2_sync(pilotRegionLen+1 : pilotRegionLen+dataLen);

pilot_rx3_slot1_time = rx3_sync(1                : pilotLen);
pilot_rx3_slot2_time = rx3_sync(pilotLen+1       : 2*pilotLen);
pilot_rx3_slot3_time = rx3_sync(2*pilotLen+1     : 3*pilotLen);
data_rx3_time        = rx3_sync(pilotRegionLen+1 : pilotRegionLen+dataLen);

% === 2) TX pilot: time->freq ===
pilot_tx_time_grid = reshape(s_vec_tx_pilot, N + Ncp, numsym_pilot);
pilot_tx_time_noCP = pilot_tx_time_grid(Ncp+1:end, :);
pilot_TX_freq_full = fft(pilot_tx_time_noCP, N, 1) / sqrt(N);
pilot_TX_freq_sub  = pilot_TX_freq_full(posLower, :);

% === 3) RX pilot: time->freq (subcarriers) ===
pilot_RX1_1 = slot_fft_sub(pilot_rx1_slot1_time, N, Ncp, numsym_pilot, posLower);
pilot_RX1_2 = slot_fft_sub(pilot_rx1_slot2_time, N, Ncp, numsym_pilot, posLower);
pilot_RX1_3 = slot_fft_sub(pilot_rx1_slot3_time, N, Ncp, numsym_pilot, posLower);

pilot_RX2_1 = slot_fft_sub(pilot_rx2_slot1_time, N, Ncp, numsym_pilot, posLower);
pilot_RX2_2 = slot_fft_sub(pilot_rx2_slot2_time, N, Ncp, numsym_pilot, posLower);
pilot_RX2_3 = slot_fft_sub(pilot_rx2_slot3_time, N, Ncp, numsym_pilot, posLower);

pilot_RX3_1 = slot_fft_sub(pilot_rx3_slot1_time, N, Ncp, numsym_pilot, posLower);
pilot_RX3_2 = slot_fft_sub(pilot_rx3_slot2_time, N, Ncp, numsym_pilot, posLower);
pilot_RX3_3 = slot_fft_sub(pilot_rx3_slot3_time, N, Ncp, numsym_pilot, posLower);

% === 4) RX data: time->freq ===
data_rx1 = slot_fft_sub(data_rx1_time, N, Ncp, numsym_data, posLower);
data_rx2 = slot_fft_sub(data_rx2_time, N, Ncp, numsym_data, posLower);
data_rx3 = slot_fft_sub(data_rx3_time, N, Ncp, numsym_data, posLower);

% === 5) H 估计（与你 RX 一致） ===
eps_val = 1e-12;

R11 = pilot_RX1_1 ./ (pilot_TX_freq_sub + eps_val);
R12 = pilot_RX1_2 ./ (pilot_TX_freq_sub + eps_val);
R13 = pilot_RX1_3 ./ (pilot_TX_freq_sub + eps_val);

R21 = pilot_RX2_1 ./ (pilot_TX_freq_sub + eps_val);
R22 = pilot_RX2_2 ./ (pilot_TX_freq_sub + eps_val);
R23 = pilot_RX2_3 ./ (pilot_TX_freq_sub + eps_val);

R31 = pilot_RX3_1 ./ (pilot_TX_freq_sub + eps_val);
R32 = pilot_RX3_2 ./ (pilot_TX_freq_sub + eps_val);
R33 = pilot_RX3_3 ./ (pilot_TX_freq_sub + eps_val);

H11 = mean(R11, 2); H12 = mean(R12, 2); H13 = mean(R13, 2);
H21 = mean(R21, 2); H22 = mean(R22, 2); H23 = mean(R23, 2);
H31 = mean(R31, 2); H32 = mean(R32, 2); H33 = mean(R33, 2);

numSC = numel(posLower);
H_all = [H11 H12 H13 H21 H22 H23 H31 H32 H33];
H = permute(reshape(H_all, numSC, 3, 3), [2 3 1]);   % 3x3xSC

% === 6) ZF 等化（逐子载波） ===
numsym = numsym_data;
Xhat = zeros(numSC, numsym, 3);

for k = 1:numSC
    Yk = [data_rx1(k,:); data_rx2(k,:); data_rx3(k,:)].';
    Xk = Yk / H(:,:,k);
    Xhat(k,:,:) = Xk;
end

eq_data1 = squeeze(Xhat(:,:,1));
eq_data2 = squeeze(Xhat(:,:,2));
eq_data3 = squeeze(Xhat(:,:,3));

% === 7) 噪声方差估计（你 RX 用的是 1e-2 floor，我建议扫 SNR 时改小一点） ===
[noiseVar1, noiseVar2, noiseVar3] = estimate_noisevar_from_pilots( ...
    H, pilot_RX1_1, pilot_RX2_1, pilot_RX3_1, ...
       pilot_RX1_2, pilot_RX2_2, pilot_RX3_2, ...
       pilot_RX1_3, pilot_RX2_3, pilot_RX3_3, ...
    pilot_TX_freq_sub);

% === 8) LLR ===
llr1 = qamdemod(eq_data1(:), M, 'gray', 'OutputType','llr', 'UnitAveragePower', false, 'NoiseVariance', noiseVar1);
llr2 = qamdemod(eq_data2(:), M, 'gray', 'OutputType','llr', 'UnitAveragePower', false, 'NoiseVariance', noiseVar2);
llr3 = qamdemod(eq_data3(:), M, 'gray', 'OutputType','llr', 'UnitAveragePower', false, 'NoiseVariance', noiseVar3);

% === 9) 联合 LLR 堆叠 + 全局解交织 + 联合 LDPC 解码（与你 RX 一致） ===
[~, cfgDec_rx] = getProtoMatrix(648, 540);

L1 = reshape(llr1, bitsPerSymbol, []);
L2 = reshape(llr2, bitsPerSymbol, []);
L3 = reshape(llr3, bitsPerSymbol, []);

L_stack = [L1; L2; L3];
llr_combined = L_stack(:);

rng(42);  % 必须与 TX 一致（你 TX 强制 rng(42)）
len_syms = length(llr_combined) / bitsPerSymbol;
perm_idx = randperm(len_syms);

llr_per_symbol = reshape(llr_combined, bitsPerSymbol, []);
inv_perm = zeros(1, len_syms);
inv_perm(perm_idx) = 1:len_syms;

llr_restored = llr_per_symbol(:, inv_perm);
llr_final_stream = llr_restored(:);

decBits_combined = decodeLdpcStreamSimple(llr_final_stream, K, blockLength, cfgDec_rx, numIter);

totalInfoBits = length(txBits1) + length(txBits2) + length(txBits3);
if length(decBits_combined) < totalInfoBits
    decBits_combined = [decBits_combined; zeros(totalInfoBits - length(decBits_combined), 1)];
end
decBits_combined = decBits_combined(1:totalInfoBits);

bits_reshaped = reshape(decBits_combined, 3, []);
decBits1 = bits_reshaped(1,:).';
decBits2 = bits_reshaped(2,:).';
decBits3 = bits_reshaped(3,:).';

% === 10) BER ===
len1 = min(length(decBits1), length(txBits1));
len2 = min(length(decBits2), length(txBits2));
len3 = min(length(decBits3), length(txBits3));

err1 = nnz(txBits1(1:len1) ~= decBits1(1:len1));
err2 = nnz(txBits2(1:len2) ~= decBits2(1:len2));
err3 = nnz(txBits3(1:len3) ~= decBits3(1:len3));

err3 = [err1 err2 err3];
tot3 = [len1 len2 len3];

end

%% ===== helper funcs (从你代码抽出来) =====
function Xsub = slot_fft_sub(slot_time_vec, N, Ncp, numsym, posLower)
    grid = reshape(slot_time_vec, N + Ncp, numsym);
    noCP = grid(Ncp+1:end,:);
    full = fft(noCP, N, 1) / sqrt(N);
    Xsub = full(posLower,:);
end

function [noiseVar1, noiseVar2, noiseVar3] = estimate_noisevar_from_pilots( ...
    H, RX1_1, RX2_1, RX3_1, RX1_2, RX2_2, RX3_2, RX1_3, RX2_3, RX3_3, TX)

    numSC = size(TX,1);
    Np    = size(TX,2);

    rxeq1 = zeros(numSC, Np);
    rxeq2 = zeros(numSC, Np);
    rxeq3 = zeros(numSC, Np);

    for ksc = 1:numSC
        Hk = H(:,:,ksc);

        Y1 = [RX1_1(ksc,:); RX2_1(ksc,:); RX3_1(ksc,:)].';
        X1 = Y1 / Hk; rxeq1(ksc,:) = X1(:,1).';

        Y2 = [RX1_2(ksc,:); RX2_2(ksc,:); RX3_2(ksc,:)].';
        X2 = Y2 / Hk; rxeq2(ksc,:) = X2(:,2).';

        Y3 = [RX1_3(ksc,:); RX2_3(ksc,:); RX3_3(ksc,:)].';
        X3 = Y3 / Hk; rxeq3(ksc,:) = X3(:,3).';
    end

    diff1 = rxeq1 - TX;
    diff2 = rxeq2 - TX;
    diff3 = rxeq3 - TX;

    noiseVar_raw = [ ...
        mean(abs(diff1(:)).^2, 'omitnan'), ...
        mean(abs(diff2(:)).^2, 'omitnan'), ...
        mean(abs(diff3(:)).^2, 'omitnan') ];

    pilotPower = mean(abs(TX(:)).^2, 'omitnan');

    % 扫 SNR 时建议更小，避免高 SNR "地板"
    noiseFloor = 1e-6 * pilotPower;

    noiseVar = max(noiseVar_raw, noiseFloor);
    noiseVar(~isfinite(noiseVar)) = noiseFloor;

    noiseVar1 = noiseVar(1);
    noiseVar2 = noiseVar(2);
    noiseVar3 = noiseVar(3);
end

function infoBitsVec = decodeLdpcStreamSimple(llrStream, K, blockLength, cfgDec, numIter)
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

% DCO-OFDM Transmitter WDM interleaving
% Author: Akatsuki Sky

clc; 
clear; 
close all;
rng(42);
% ----------------- 路径处理 -----------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

% ----------------- 基本参数 -----------------
numChannels      = 3;    
M                = 4;   
bitsPerSymbol    = log2(M);

N                = 256; 
Ncp              = 32;   
numsym_data      = 200;  % 数据 OFDM 符号个数
numsym_pilot     = 100;  % 导频 OFDM 符号个数

numDataSubcarriers   = N/2 - 1;            
symbolsPerTX_needed  = numDataSubcarriers * numsym_data;
bitsPerTX            = symbolsPerTX_needed * bitsPerSymbol;
numBitsPerTX         = bitsPerTX;
totalQAMsymbolsNeeded= numChannels * symbolsPerTX_needed;

% 剪切门限（DCO-OFDM）
mu_bottom_dB = 10;
mu_top_dB    = 10;
mu_bottom    = -(10^(mu_bottom_dB/20))*sqrt(2);
mu_top       =  (10^(mu_top_dB/20))*sqrt(2);

posLower = 2:(N/2);    
posUpper = (N/2+2):N;  

% LDPC 参数
[cfgEnc, cfgDec] = getProtoMatrix(648, 540);
K           = cfgEnc.NumInformationBits;
blockLength = cfgDec.BlockLength;
numIter     = 25;      

% =====================================================
% 1. 三路比特生成 + 联合交织 + 联合 LDPC 编码 + 符号分发
% =====================================================

totalPhysicalBits = totalQAMsymbolsNeeded * bitsPerSymbol; 
numBlocksFit = floor(totalPhysicalBits / blockLength);
totalInfoBits_System = numBlocksFit * K;
numBitsPerTX = floor(totalInfoBits_System / 3);

% 生成原始比特
txBits1 = randi([0 1], numBitsPerTX, 1);
txBits2 = randi([0 1], numBitsPerTX, 1);
txBits3 = randi([0 1], numBitsPerTX, 1);

% 1. 比特级交织
bitInterleaved = reshape([txBits1.'; txBits2.'; txBits3.'], [], 1);

% 2. 联合 LDPC 编码
codedBits_combined = encodeLdpcStreamSimple(bitInterleaved, K, cfgEnc);

% 3. 填充 (Padding) 与 QAM 映射
% 编码后的比特数可能略小于物理层容量，需要补零填充，否则 reshape 会报错
numCoded = length(codedBits_combined);
if numCoded < totalPhysicalBits
    paddingLen = totalPhysicalBits - numCoded;
    codedBits_padded = [codedBits_combined; zeros(paddingLen, 1)];
else
    % 如果超出(极少情况)，则截断
    codedBits_padded = codedBits_combined(1:totalPhysicalBits);
end

qamStream_combined = qammod(codedBits_padded, M, 'gray', 'InputType', 'bit');


% 插入全局随机交织 (Interleaving)
rng(42); % 1. 强制重置种子 (必须与接收端完全一致！)
perm_idx = randperm(length(qamStream_combined));
qamStream_interleaved = qamStream_combined(perm_idx);

% 4. 符号级分发 
% 【注意】这里的输入变量必须改为 qamStream_interleaved
qam_reshaped = reshape(qamStream_interleaved, 3, []).'; 


dataSym1 = reshape(qam_reshaped(:, 1), numDataSubcarriers, numsym_data);
dataSym2 = reshape(qam_reshaped(:, 2), numDataSubcarriers, numsym_data);
dataSym3 = reshape(qam_reshaped(:, 3), numDataSubcarriers, numsym_data);
% =====================================================
% 2. 构造三路频域 OFDM 符号（Hermitian 对称） → IFFT → 加 CP → 裁剪
% =====================================================

% ------------- 信道 1 -------------
X1              = zeros(N, numsym_data);
X1(posLower,:)  = dataSym1;
X1(posUpper,:)  = flipud(conj(X1(posLower,:))); 
X1(1,:)         = 0;                             
X1(N/2+1,:)     = 0;                              

timeDomain1     = ifft(X1, N, 1) / (1/sqrt(N));  
cp1             = timeDomain1(end-Ncp+1:end, :);
s1              = [cp1; timeDomain1];      
s1(s1<mu_bottom)= mu_bottom;              
s1(s1>mu_top)   = mu_top;
s_vec_tx_data1  = s1(:);                     

% ------------- 信道 2 -------------
X2              = zeros(N, numsym_data);
X2(posLower,:)  = dataSym2;
X2(posUpper,:)  = flipud(conj(X2(posLower,:)));
X2(1,:)         = 0;
X2(N/2+1,:)     = 0;

timeDomain2     = ifft(X2, N, 1) / (1/sqrt(N));
cp2             = timeDomain2(end-Ncp+1:end, :);
s2              = [cp2; timeDomain2];
s2(s2<mu_bottom)= mu_bottom;
s2(s2>mu_top)   = mu_top;
s_vec_tx_data2  = s2(:);

% ------------- 信道 3 -------------
X3              = zeros(N, numsym_data);
X3(posLower,:)  = dataSym3;
X3(posUpper,:)  = flipud(conj(X3(posLower,:)));
X3(1,:)         = 0;
X3(N/2+1,:)     = 0;

timeDomain3     = ifft(X3, N, 1) / (1/sqrt(N));
cp3             = timeDomain3(end-Ncp+1:end, :);
s3              = [cp3; timeDomain3];
s3(s3<mu_bottom)= mu_bottom;
s3(s3>mu_top)   = mu_top;
s_vec_tx_data3  = s3(:);

% =====================================================
% 3. 导频生成
% =====================================================

pilotBits   = randi([0 1], numDataSubcarriers*numsym_pilot*bitsPerSymbol, 1);
pilotQAM    = qammod(pilotBits, M, 'gray', 'InputType', 'bit');
pilotData   = reshape(pilotQAM, numDataSubcarriers, numsym_pilot);

Xpilot              = zeros(N, numsym_pilot);
Xpilot(posLower,:)  = pilotData;
Xpilot(posUpper,:)  = flipud(conj(Xpilot(posLower,:)));
Xpilot(1,:)         = 0;
Xpilot(N/2+1,:)     = 0;

Xpilot              = Xpilot * Norcoeffi(M);

xpilot              = ifft(Xpilot, N, 1) / (1/sqrt(N));
cp_pilot            = xpilot(end-Ncp+1:end, :);
spilot              = [cp_pilot; xpilot];

spilot(spilot<mu_bottom) = mu_bottom;
spilot(spilot>mu_top)    = mu_top;

s_vec_tx_pilot      = spilot(:);

% =====================================================
% 4. 三路信道的导频时隙错时复用（no cross interleaving）
% =====================================================

pilotSlot       = length(s_vec_tx_pilot);
zeroPilot       = zeros(pilotSlot,1);


%zero = zeros(1,19600)';

% 信道 1：导频在第 1 段
tx1 = [s_vec_tx_pilot; zeroPilot;      zeroPilot;      s_vec_tx_data1];

% 信道 2：导频在第 2 段
tx2 = [zeroPilot;      s_vec_tx_pilot; zeroPilot;      s_vec_tx_data2];

% 信道 3：导频在第 3 段
tx3 = [zeroPilot;      zeroPilot;      s_vec_tx_pilot; s_vec_tx_data3];

% =====================================================
% 5. 保存参数 + 波形
% =====================================================

finalSerialSignal     = tx1;
finalSerialSignal_ch2 = tx2;
finalSerialSignal_ch3 = tx3;
upsampleRate          = 1;

pilotLen = (N + Ncp) * numsym_pilot;
dataLen  = (N + Ncp) * numsym_data;

allPath = fullfile(scriptDir, 'Knowledge.mat');

save(allPath, ...
    'numChannels','M','bitsPerSymbol','N','Ncp','numsym_pilot','numsym_data', ...
    'numDataSubcarriers','numBitsPerTX','symbolsPerTX_needed', ...
    'totalQAMsymbolsNeeded','pilotLen','dataLen', ...
    'mu_bottom_dB','mu_top_dB','mu_bottom','mu_top', ...
    'posLower','K','blockLength','numIter', ...
    'txBits1','txBits2','txBits3', ...
    's_vec_tx_data1','s_vec_tx_data2','s_vec_tx_data3', ...
    's_vec_tx_pilot','Xpilot', ...
    'finalSerialSignal','finalSerialSignal_ch2','finalSerialSignal_ch3', ...
    'upsampleRate');

% =====================================================
% 6. 导出 CSV（每一路一个文件）
% =====================================================

csvwrite(fullfile(scriptDir, 'tx1_simple.csv'), tx1);
csvwrite(fullfile(scriptDir, 'tx2_simple.csv'), tx2);
csvwrite(fullfile(scriptDir, 'tx3_simple.csv'), tx3);

function codedBits = encodeLdpcStreamSimple(bitStream, K, cfgEnc)
    numBlocks = ceil(numel(bitStream) / K);
    padLen    = numBlocks*K - numel(bitStream);
    if padLen > 0
        bitStream = [bitStream; zeros(padLen,1)];
    end
    infoMatrix = reshape(bitStream, K, numBlocks);
    codeword   = ldpcEncode(infoMatrix, cfgEnc, "OutputFormat", "whole");
    codedBits  = codeword(:);              
end

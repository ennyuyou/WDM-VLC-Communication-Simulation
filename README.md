# WDM-VLC-Communication-Simulation

MATLAB simulation for a multi-channel Wavelength Division Multiplexing Visible Light Communication (WDM-VLC) system utilizing advanced signal processing and error correction.

---

## 概要 / Project Overview

### 日本語
本プロジェクトは、LEDと蛍光体ファイバアンテナを用いたマルチチャネル波長分割複用可視光通信（WDM-VLC）システムのMATLABシミュレーションである。
低密度パリティ検査（LDPC）符号および離散マルチトーン（DCO-OFDM）変調を統合し、信号処理技術の最適化を行うことで、高速かつ低誤り率のデータ伝送を達成している。
サンプリングレート 25 MSa/s において、総スループット 112.5 Mbps を達成し、ビット誤り率（BER）を硬判決前向誤り訂正（HD-FEC）限度以下に抑えるシステム検証を可能とする。

### English
This project implements a MATLAB simulation platform for a multi-channel Wavelength Division Multiplexing Visible Light Communication (WDM-VLC) system using LEDs and fluorescent fiber antennas. 
By integrating Low-Density Parity-Check (LDPC) coding and Direct-Current-Biased Optical OFDM (DCO-OFDM), the simulation optimizes data transmission efficiency. At a sampling rate of 25 MSa/s, the system achieves a total throughput of 112.5 Mbps while maintaining the Bit Error Rate (BER) below the hard-decision forward error correction (HD-FEC) threshold.

---

## 主な機能 / Key Features

* **Error Correction (誤り訂正):** Advanced LDPC encoding and decoding configurations to combat channel noise.
* **Modulation (変調方式):** DCO-OFDM modulation with dynamic bit-level/symbol-level interleaving.
* **WDM Channel Simulation (多チャネルシミュレーション):** Simulates a 3-channel optical communication environment, factoring in signal attenuation and fluorescent response.
* **Performance Analysis (性能評価):** Evaluation of system performance via Bit Error Rate (BER) curves and constellation diagrams.

---

## 開発環境 / Environment
* MATLAB (R2022a or later recommended)
* Communications Toolbox

---

## 成果物展示 / Results

> 💡 **Tip for Interviewers:** > 詳しいシミュレーション結果（誤り率曲線やコンスタレーション図）は、ソースコードを実行することで確認できます。また、研究の詳細はJAISTでの修士論文に基づいています。
> (Detailed simulation results, including BER curves, can be reproduced by running the main scripts. The logic is based on my Master's thesis at JAIST.)

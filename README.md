# Viterbi Decoder for 4-PAM over ISI Channels

## Overview
This project implements a **Viterbi decoder** for symbol detection in a **4-PAM (Pulse Amplitude Modulation)** system affected by **Inter-Symbol Interference (ISI)**. The decoder is based on the **Maximum Likelihood Sequence Estimation (MLSE)** algorithm, which uses the **Viterbi Algorithm** to decode transmitted symbols in the presence of ISI.  

## Theoretical Background

###  **Channel Model**
The received signal at time `k` can be modeled as:  
\$y_k = \sum_{i=0}^{\mu} h_i x_{k-i} + n_k$

where:  
- $( x_k $) is the transmitted 4-PAM symbol.  
- $( h_i $) represents the ISI channel coefficients.  
- $( \mu $) is the **channel memory** (ISI order).  
- $( n_k $) is additive white Gaussian noise (AWGN) with variance $( \sigma^2 $).  

### **Viterbi Algorithm for ISI Equalization**
The **Viterbi decoder** is used for MLSE by:  
1. **State Representation**: Each state represents a possible sequence of past transmitted symbols.  
2. **Trellis Construction**: Transition probabilities are computed using the channel model.  
3. **Path Metric Computation**: The decoder selects the path with the **minimum accumulated Euclidean distance**.  
4. **Traceback & Decision**: The most likely sequence is estimated based on the **survivor path**.

## Features
✔️ Implements **Viterbi decoding** for **4-PAM ISI equalization**.  
✔️ Supports **different ISI channel memory depths (μ)**.  
✔️ Computes **Symbol Error Probability (SEP) vs. SNR curves**.  
✔️ **Optimized traceback depth** for efficient sequence estimation.  
✔️ **Trellis visualization** for debugging & analysis.  


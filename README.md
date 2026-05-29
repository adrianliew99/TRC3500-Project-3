# TRC3500-Project-3
Project 3 github repositories for students from TRC3500 Wed 9am Group A4. STM32-based breath-rate estimation system using conductive-rubber strain sensing, TMP61 thermistor airflow sensing, digital signal processing, and adaptive sensor-fusion for robust respiratory monitoring.

# Team Member Introduction
### Member 1: The claude cooker 
Lim Yi Woeh

### Member 2: The night shifter 
Loke Wai Hon 
Who work start from 2am til morning 8am he's crazy

### Member 3: The jit gor
Lee Zong Xuan 

### Member 4: The diplomatic officer AKA External Lialson officer 
##### no lah just sometimes go steal component and source component from a lot friend friend

Liew Chong Jun

### Member 5: The emotional value provider
Benjamin Ooi Jing Yew 

# Thanos
<img width="360" height="364" alt="tanwenshandrinkingcoffee" src="https://github.com/user-attachments/assets/ba95a96d-7a70-47f2-b161-804d2dd7d203" />


## Project Description

This repository contains the source code, circuit design files, data-processing scripts, and documentation for **TRC3500 Project 3: Breath Rate Estimation**. The project develops a multi-sensor respiratory monitoring prototype that estimates human breathing rate using two complementary sensing methods.

The system combines a **conductive-rubber strain sensor** and a **TMP61 thermistor**. The conductive rubber sensor is mounted around the lower sternum to detect chest expansion during breathing, while the thermistor is positioned near the nostril to detect temperature changes between inhaled and exhaled airflow. Both sensor outputs are conditioned using op-amp interface circuits and acquired by an **STM32L432KC** through ADC and DMA.

The collected data is streamed to a Python-based signal-processing pipeline, where baseline drift removal, smoothing, Hilbert-phase breath detection, and event-level sensor fusion are applied. The fusion algorithm adaptively selects the more reliable sensor depending on the breathing condition, producing a unified breath-rate estimate that is more robust than relying on a single sensor alone.

The project includes experimental evaluation using metronome-guided breathing as ground truth, with performance assessed using breath-rate error, RMSE, error histograms, and just-noticeable-difference analysis.


# Appreciation and Credit
### Jensen Huang 
<img width="225" height="225" alt="jensenhuanggreeneye" src="https://github.com/user-attachments/assets/835328e2-3956-4ba5-a1f4-75cb97a72e29" />

### YiLong Ma
<img width="2160" height="3840" alt="yilongma" src="https://github.com/user-attachments/assets/69504033-3c72-41f1-bb9c-7b100cf6264c" />

### Zhong Xina
<img width="1700" height="2328" alt="ZhongXina" src="https://github.com/user-attachments/assets/8766107d-808c-43f3-883e-dbb139d9ecbf" />


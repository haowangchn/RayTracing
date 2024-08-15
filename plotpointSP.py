import numpy as np
import matplotlib.pyplot as plt

# 参数定义
c = 1500  # 声速，单位：m/s
A = 1     # 声源振幅

# 定义声源和镜像源的位置
r_S = 100        # 声源到接收点的距离
r_S_prime = 116.6 # 镜像源到接收点的距离

# 频率范围
f = np.linspace(10, 1000, 500)  # 频率范围，从10 Hz到1000 Hz

# 计算波数和角频率
k = (2 * np.pi * f) / c
omega = 2 * np.pi * f

# 计算声压幅度
# 声压的幅度公式
amplitude_S = A / r_S
amplitude_S_prime = A / r_S_prime
amplitude = amplitude_S - amplitude_S_prime

# 由于 amplitude 是常数，所以计算其绝对值时需要对每个频率的值进行计算
amplitude = np.abs(amplitude_S - amplitude_S_prime)

# 绘制图像
plt.figure(figsize=(8, 6))
plt.plot(f, amplitude * np.ones_like(f), label='Sound Pressure Amplitude')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Sound Pressure Amplitude vs. Frequency')
plt.grid(False)
plt.legend()
plt.show()

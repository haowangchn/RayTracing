import numpy as np
import matplotlib.pyplot as plt

# 参数定义
c = 1500  # 声速，单位：m/s
A = 1     # 声源振幅
r_S = 100        # 声源到接收点的距离
r_S_prime = 116.6 # 镜像源到接收点的距离

# 频率范围
f = np.linspace(10, 1000, 50)  # 频率范围，从10 Hz到1000 Hz
omega = 2 * np.pi * f  # 角频率

# 计算波数 k
k = omega / c

# 计算声压 (复数形式)
p_S = A / r_S * np.exp(-1j * k * r_S)  # 计算声源的声压
p_S_prime = A / r_S_prime * np.exp(-1j * k * r_S_prime)  # 计算镜像源的声压

# 相干声压总和
p_total = p_S + p_S_prime

# 计算声压幅值
amplitude = np.abs(p_total)

# 绘制不带光滑拟合的折线图
plt.figure(figsize=(10, 6))
plt.plot(f, amplitude, '-', label='Coherent Sound Pressure Amplitude')  # 使用 '-' 表示折线
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('Coherent Sound Pressure Amplitude vs. Frequency')
plt.grid(False)
plt.legend()

# 保存图像到文件
plt.savefig('coherent_sound_pressure_amplitude.eps')
# plt.show()  # 显示图像

# 将频率和声压幅值保存到文件中
data = np.column_stack((f, amplitude))
np.savetxt('data1.txt', data, header='Frequency (Hz)  Amplitude', fmt='%10.4f', delimiter='\t')

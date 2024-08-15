import numpy as np
import matplotlib.pyplot as plt

# 读取数据
r = np.loadtxt('output_r.txt')
z = np.loadtxt('output_z.txt')
psi = np.loadtxt('output_psi.txt')

# 将psi数据转换为复数
psi = psi[:, 0] + 1j * psi[:, 1]
psi = psi.reshape(len(z), len(r))

# 选择部分数据进行绘图
takez = slice(0, len(z), 5)
taker = slice(1, len(r), 10)

zt = z[takez]
rt = r[taker]
psit = psi[takez, taker]

# 计算 Hankel 函数
k0 = 2 * np.pi * 50 / 1500
hank = np.sqrt(2 / (np.pi * k0)) * np.exp(1j * (k0 * rt - np.pi / 4)) / np.sqrt(rt)

# 计算传输损失
tl = 20 * np.log10(np.abs(psit * hank))

# 绘图
plt.pcolor(rt, zt, tl, shading='auto')
plt.colorbar()
plt.clim(-100, -60)
plt.xlabel('Range (m)')
plt.ylabel('Depth (m)')
plt.title('PE intensity')
plt.show()

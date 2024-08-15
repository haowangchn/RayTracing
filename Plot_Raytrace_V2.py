# test

import matplotlib.pyplot as plt
import numpy as np

# 从文件读取数据
data = np.loadtxt('test_RayTrace.txt', skiprows=2)  # 跳过第一行标题

# 提取数据列
depth = data[:, 0]
horizontal_length = data[:, 1]
sound_speed = data[:,2]
A0 = data[:,5]

# 设置图形大小和布局
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6))

# 第一个子图：水深和水平长度
ax1.plot(horizontal_length, depth, linestyle='-', color='b', linewidth=1)
ax1.set_title('Ray Tracing Results')
ax1.set_xlabel('Horizontal Length')
ax1.set_ylabel('Depth')
ax1.invert_yaxis()  # 翻转Y轴，箭头朝下
# ax1.invert_xaxis()

# 第二个子图：水深和声速
ax2.plot(sound_speed, depth, linestyle='-', color='g',linewidth=1)
ax2.set_title('Sound Speed Profile')
ax2.set_xlabel('Sound Speed')
ax2.set_ylabel('Depth')
ax2.invert_yaxis()  # 翻转Y轴，箭头朝下

# 第三个子图：幅值--水深
ax3.plot(A0[:50], depth[:50], linestyle='-', color='r', linewidth=1)
ax3.set_title('Amplitude')
ax3.set_xlabel('A0')
ax3.set_ylabel('Depth')
ax3.invert_yaxis()  # 翻转Y轴，箭头朝下

# 显示图形
plt.tight_layout()
plt.show()

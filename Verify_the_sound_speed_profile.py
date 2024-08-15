# verify the sound speed profile
# 验证结果：符合算例

import matplotlib.pyplot as plt
import numpy as np

# 从文件读取数据
data = np.loadtxt('test_RayTrace.txt', skiprows=2)  # 跳过第一行标题
data_verify = np.loadtxt('sound_speed_profile.txt', skiprows = 1)

# 提取数据列
depth = data[:, 0]
sound_speed = data[:,2]

# 用以验证的标准数据集
depth_verify = data_verify[:,0]
sound_speed_verify = data_verify[:,1]

# 设置图形大小和布局
fig, (ax1) = plt.subplots(1, 1, figsize=(5, 6))

# 声速剖面
ax1.plot(sound_speed_verify, depth_verify, linestyle='-', color='r', linewidth=1, label='Standard Profile', markersize=2)
ax1.plot(sound_speed, depth, linestyle='-', color='g', linewidth=4, label='Calculated Profile', markersize=2, alpha=0.3)

ax1.set_title('Sound Speed Profile')
ax1.set_xlabel('Sound Speed')
ax1.set_ylabel('Depth')
ax1.invert_yaxis()  # 翻转Y轴，箭头朝下


# 显示图形
plt.tight_layout()
plt.show()

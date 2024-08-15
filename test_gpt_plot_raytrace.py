import numpy as np
import matplotlib.pyplot as plt

# 读取数据文件
data = np.loadtxt('test_RayTrace.txt', skiprows=1)

# 提取数据
ray_number = data[:, 0].astype(int)   # 声线编号
iteration = data[:, 1].astype(int)    # 迭代次数
z = data[:, 2]                        # 深度数据
r = data[:, 3]                        # 水平距离数据
c = data[:, 4]                        # 声速数据

# 获取唯一的声线编号
unique_rays = np.unique(ray_number)

# 创建画布
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# 绘制声线图
for ray in unique_rays[:2]:  # 只绘制前两条声线
    mask = (ray_number == ray)
    ax1.plot(r[mask], z[mask], label=f'Ray {ray}')
ax1.set_xlabel('Horizontal Distance (r)')
ax1.set_ylabel('Depth (z)')
ax1.invert_yaxis()  # 将深度翻转，使得0在最上方
ax1.set_title('Ray Paths')
ax1.legend()

# 绘制声速剖面图
ax2.plot(c, z, marker='o', linestyle='-', color='b')
ax2.set_xlabel('Sound Speed (c)')
ax2.set_ylabel('Depth (m)')
ax2.invert_yaxis()  # 将深度翻转，使得0在最上方
ax2.set_title('Sound Speed Profile')

# 调整布局，显示图形
plt.tight_layout()
plt.show()

import pandas as pd
import matplotlib.pyplot as plt

# 读取文件
data = pd.read_csv('SoundProfile.txt', skiprows=2, delimiter='\s+', names=['Depth', 'SoundSpeed'])

# 绘制图形
plt.figure(figsize=(5, 6))
plt.plot(data['SoundSpeed'], data['Depth'], marker='o', linestyle='-')

plt.xlabel('Sound Speed (m/s)')
plt.ylabel('Depth (m)')
plt.title('Sound Profile')

plt.gca().invert_yaxis()  # 反转 y 轴，使深度从上到下增加
plt.grid(False)
plt.show()

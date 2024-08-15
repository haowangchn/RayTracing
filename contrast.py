import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager

# 设置中文字体，确保中文能够正常显示
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体字体
plt.rcParams['axes.unicode_minus'] = False    # 正常显示负号

# 读取数据文件，跳过第一行
data1 = pd.read_csv('data1.txt', delimiter='\t', skiprows=1, header=None)
data2 = pd.read_csv('data2.txt', delimiter='\t', skiprows=1, header=None)

# 手动指定列名
data1.columns = ['Frequency (Hz)', 'Sound Pressure Magnitude (Pa)']
data2.columns = ['Frequency (Hz)', 'Sound Pressure Magnitude (Pa)']

# 绘制图表
plt.figure(figsize=(10, 6))

# 绘制data1
plt.plot(data1['Frequency (Hz)'], data1['Sound Pressure Magnitude (Pa)'], label='数值解', marker='o')

# 绘制data2
plt.plot(data2['Frequency (Hz)'], data2['Sound Pressure Magnitude (Pa)'], label='Bellhop解', marker='s')

# 添加图例、标题和标签
plt.legend()
plt.title('Frequency vs Sound Pressure Magnitude')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Sound Pressure Magnitude (Pa)')

# 显示图表
plt.show()

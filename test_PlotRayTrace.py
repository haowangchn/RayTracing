import matplotlib.pyplot as plt

# 读取声线轨迹文件
with open('test_RayTrace.txt', 'r') as file:
    lines = file.readlines()

depth_ray_trace = []
horizontal_length = []
current_bundle = []

for line in lines:
    if line.startswith('Ray Bundle'):
        if current_bundle:
            depth_ray_trace.append(current_bundle)
            current_bundle = []
    else:
        data = line.split()
        current_bundle.append((float(data[0]), float(data[1])))

if current_bundle:
    depth_ray_trace.append(current_bundle)

# 绘制图形
fig, ax = plt.subplots(figsize=(10, 6))

# 绘制每条声线束的轨迹
for bundle in depth_ray_trace:
    depths, lengths = zip(*bundle)
    ax.plot(lengths, depths)

# 标注点 (rfp, zfp)
rfp = 1.0
zfp = 1000.0
ax.plot(rfp, zfp, 'ro', label='Field Point (rfp, zfp)')  # 使用红色圆点标注
ax.annotate(f'({rfp}, {zfp})', (rfp, zfp), textcoords="offset points", xytext=(0,10), ha='center')  # 添加标注文字

ax.invert_yaxis()  # 纵坐标从上到下
ax.xaxis.tick_top()  # 横轴标注在上框
ax.xaxis.set_label_position('top')  # 横轴标签在上框
ax.set_xlabel('Horizontal Length')
ax.set_ylabel('Depth')
ax.set_title('Ray Trace vs Depth')
ax.legend()  # 显示图例

plt.tight_layout()
plt.show()

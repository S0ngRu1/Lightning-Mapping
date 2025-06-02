import math

def angle_with_y_axis(x1, y1, x2, y2):
    """计算从点(x1, y1)到点(x2, y2)的向量与Y轴正方向的夹角"""
    # 计算向量
    vector = (x2 - x1, y2 - y1)
    # 计算向量的模长
    magnitude = math.sqrt(vector[0]**2 + vector[1]**2)
    # 计算向量与Y轴正方向单位向量的点积
    dot_product = vector[1]
    # 计算余弦值
    cos_theta = dot_product / magnitude
    # 计算角度
    angle = math.acos(cos_theta)
    # 将弧度转换为度
    angle_degrees = math.degrees(angle)
    return angle_degrees


# 假设点A, B, C的坐标分别为(ax, ay), (bx, by), (cx, cy)
ax, ay = 0, 0
bx, by = -2.0622,41.5985
cx, cy = 28.4316,23.5237

# 计算每条边与Y轴的夹角
angle_AB = angle_with_y_axis(ax, ay, bx, by)
angle_BC = angle_with_y_axis(bx, by, cx, cy)
angle_CA = angle_with_y_axis(cx, cy, ax, ay)

print(f"Angle between AB and Y-axis: {angle_AB} degrees")
print(f"Angle between BC and Y-axis: {angle_BC} degrees")
print(f"Angle between CA and Y-axis: {angle_CA} degrees")
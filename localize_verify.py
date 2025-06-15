# -*- coding: utf-8 -*-
# @Time : 2025/4/27 22:57
# @Author : CSR
# @File : localize_verify.py

import numpy as np
import math

# 常量
c = 299792458  # 光速，单位：米/秒
eps = 1e-10  # 用于避免除以零的小量
rad2deg = 180 / np.pi
deg2rad = np.pi / 180

# 站点位置 (假设在同一个局部笛卡尔坐标系下)
# 假设 +X 为东, +Y 为北, +Z 为上
yld_sit = np.array([0, 0, 0], dtype=float)
chj_sit = np.array([1991, -7841.2, 0], dtype=float)  # yld 相对于 chj 的位置


# ----------------------------------------------------------------------
# 辅助函数：坐标系转换

def direction_vector_from_angles(azimuth_from_N_deg, elevation_from_horz_deg):
    """
    将标准方位角（从北+Y向东+X）和仰角（从水平XY平面向上）转换为
    单位方向向量（+X东, +Y北, +Z上）。

    这是根据你的代码中 sph2cart 的用法推断的转换方式：
    sph2cart(deg2rad(90-az), deg2rad(el), 1)
    假设 sph2cart(az_rad, el_rad, r) 返回 [r*cos(el)*cos(az), r*cos(el)*sin(az), r*sin(el)]
    其中 el 是仰角从XY平面向上，az 是从+X向+Y。
    那么你的调用使用了 az_rad = deg2rad(90-azimuth_from_N_deg) 作为从X的角度，
    el_rad = deg2rad(elevation_from_horz_deg) 作为从XY平面的仰角。
    这确实是标准的方向向量转换到 +X东, +Y北, +Z上 坐标系的方式。
    """
    az_rad_for_sph = deg2rad(90 - azimuth_from_N_deg)  # 从+X向+Y的角度
    el_rad_for_sph = deg2rad(elevation_from_horz_deg)  # 从XY平面向上仰角

    # 模拟 sph2cart 的行为
    x = np.cos(el_rad_for_sph) * np.cos(az_rad_for_sph)
    y = np.cos(el_rad_for_sph) * np.sin(az_rad_for_sph)
    z = np.sin(el_rad_for_sph)
    return np.array([x, y, z])


def cart2sph_standard(vector):
    """
    将笛卡尔向量 [x, y, z] (+X东, +Y北, +Z上) 转换为
    标准方位角（度，从+Y北向+X东）和仰角（度，从XY平面向上）。
    """
    x, y, z = vector
    horizontal_distance = np.sqrt(x ** 2 + y ** 2)

    # 计算仰角（从水平面向上）
    elevation_rad = np.arctan2(z, horizontal_distance)
    elevation_deg = elevation_rad * rad2deg

    # 计算方位角（从+Y向+X）
    azimuth_rad = np.arctan2(x, y)  # arctan2(y, x) in typical math libs for angle from +X
    # but often arctan2(x, y) for angle from +Y in navigation
    azimuth_deg = azimuth_rad * rad2deg
    # 确保方位角在 0-360 度范围内
    azimuth_deg = (azimuth_deg + 360) % 360

    return azimuth_deg, elevation_deg


# ----------------------------------------------------------------------
# 模拟数据生成 (前向仿真第一部分)

def simulate_observation(S_true, station_sit):
    """
    模拟从站点 station_sit 观测到 S_true 位置闪电的角度。
    返回标准方位角（从北向东）和仰角（从水平面）。
    """
    vector_to_target = S_true - station_sit
    azimuth, elevation = cart2sph_standard(vector_to_target)
    return azimuth, elevation


# ----------------------------------------------------------------------
# 实现你的代码的核心定位逻辑 (模拟你的函数)

def calculate_sub_S_user_logic(yld_sit, chj_sit, yld_azimuth, yld_elevation, chj_azimuth, chj_elevation):
    """
    根据你的代码片段实现核心的 sub_S 计算逻辑。
    注意：这里的输入 azimuth 和 elevation 应该是你的代码所期望的格式
    （我们假设是标准方位角/仰角，然后通过 direction_vector_from_angles 处理）。
    """
    p = chj_sit - yld_sit

    # 根据你的代码，使用 direction_vector_from_angles 来获取 A1 和 A2
    A1 = direction_vector_from_angles(yld_azimuth, yld_elevation)
    A2 = direction_vector_from_angles(chj_azimuth, chj_elevation)

    C = np.cross(A1, A2)
    norm_C = np.linalg.norm(C)

    if norm_C < eps:
        # Rays are parallel or anti-parallel, can't form a plane for C
        # In a real system, you might handle this (e.g., return None or signal error)
        print("警告: 射线平行或反平行，无法计算唯一解。")
        return None

    c_unit = C / norm_C

    # 构建矩阵 M
    # M = [A1', -A2', c_unit']  -- 这里是列向量组成的矩阵
    M = np.vstack((A1, -A2, c_unit)).T  # Vstack builds rows, then transpose to get columns

    # 使用克莱姆法则求解 M * [R1_value, R2_value, R3_value]' = p'
    # 等价于解线性系统 M @ [R1_value, R2_value, R3_value]' = p
    try:
        # 使用 numpy 的线性方程求解器更稳定和高效
        coeffs = np.linalg.solve(M, p)
        R1_value, R2_value, R3_value_raw = coeffs
    except np.linalg.LinAlgError:
        print("警告: 矩阵不可逆，无法求解线性系统。")
        return None

    R1 = R1_value * A1
    R2 = R2_value * A2
    # 根据你的代码 R3 的计算方式，注意这里使用了 R3_value_raw
    # R3 = R3_value_raw / norm_C * C # This is equivalent to R3_value_raw * c_unit
    R3 = R3_value_raw * c_unit  # More direct implementation of R3_value_raw * c_unit

    # *** 注意：你的原始代码在计算 R3 后有一行 R3 = -R3; ***
    # 这个操作与求解的线性系统 R1*A1 - R2*A2 + R3*c_unit = p 不符。
    # 如果要严格模拟你的代码，应该包含 R3 = -R3。
    # 但从几何意义上看，这可能是一个错误。
    # 我们先按照求解的系统计算 R3，不进行额外的取反。
    # 如果需要严格模拟你的原始代码，请取消下一行的注释：
    # R3 = -R3 # <= 模拟原始代码中的 R3 取反

    sub_S = None
    # 根据你的代码，使用不同的公式计算 sub_S
    # 注意公式中使用了 R1_value, R2_value, R3 向量, p 向量
    # 这些公式的几何意义不明确，是需要验证的部分
    if R1_value + R2_value == 0:  # Avoid division by zero if both are zero
        # This is an edge case, perhaps return None or handle differently
        print("警告: R1_value + R2_value 为零。")
        return None

    if R1_value <= R2_value:
        # 使用第一个公式
        # sub_S = R1 + (R1_value / R2_value) * (R1_value / (R1_value + R2_value)) * R3
        # R1 是从 yld_sit 出发的向量，所以 sub_S 应该是 yld_sit + 这个向量
        sub_S = yld_sit + R1 + (R1_value / R2_value) * (R1_value / (R1_value + R2_value)) * R3
    else:
        # 使用第二个公式
        # sub_S = R2 - (R2_value / R1_value) * (R2_value / (R1_value + R1_value)) * R3 + p
        # R2 是从 chj_sit 出发的向量，p = chj_sit - yld_sit
        # S = chj_sit + R2_vector
        # 第二个公式计算的是 R2_vector - (..) * R3 + p, 这个向量再加到哪个点上？
        # 公式看起来像是在 chj_sit + R2 的基础上进行修正，或者是在 yld_sit + R1 的基础上修正
        # 假设公式意图是 S = yld_sit + (R1 vector or something derived from it)
        # 或者 S = chj_sit + (R2 vector or something derived from it)
        # 第二个公式 sub_S = R2 - (...) * R3 + p  看起来是向量形式，不是直接的点。
        # 如果 sub_S 应该是一个点，那么它应该是站点的向量 + 某个位移向量
        # 重新看你的原始代码逻辑: if else 里面直接计算 sub_S = ...
        # 没有加站点位置。这部分需要明确 sub_S 代表的是点还是位移。
        # 根据后面计算 t_chj = sqrt(sum((sub_S - chj_sit).^2))/c;
        # 这表明 sub_S 是一个表示空间点的向量。
        # 那么公式应该是在某个站点位置上加上计算出的位移。
        # 第一个公式 sub_S = R1 + ... 可能是指 S = yld_sit + (R1 + ...)
        # 第二个公式 sub_S = R2 - ... + p 可能是指 S = chj_sit + (R2 - ...) + p
        # 由于 p = chj_sit - yld_sit, chj_sit + p = chj_sit + chj_sit - yld_sit = 2*chj_sit - yld_sit (likely incorrect)
        # 让我们假设 sub_S 就是计算出的点的坐标向量，并且公式的目的是直接给出这个点。
        # 那么公式 S = R1_vec + ...; S = R2_vec - ... + p 是怎么来的？
        # R1 = R1_value * A1 是从 yld 出发的方向向量乘以距离 R1_value。
        # R2 = R2_value * A2 是从 chj 出发的方向向量乘以距离 R2_value。
        # 它们本身就表示相对于站点的位移。
        # 标准方法中 S = yld_sit + R1_value*A1 或 S = chj_sit + R2_value*A2
        # 你的公式非常规。让我们严格按照你写的向量运算来写，不猜测其几何意义，
        # 仅模拟计算过程。
        sub_S_vec = R2 - (R2_value / R1_value) * (R2_value / (R1_value + R2_value)) * R3 + p
        # 原始代码似乎没有明确将这个向量加到某个站点位置上。
        # 但后面的时间计算 dlta_t = abs(t_yld-t_chj) 使用了 sqrt(sum((sub_S - station_sit).^2))/c;
        # 这明确表明 sub_S 必须是点坐标。
        # 因此，假设你的公式计算出的 sub_S 就是空间点坐标本身。
        # 也就是说，第一个公式计算的是点坐标，第二个公式计算的也是点坐标。
        # 这意味着 R1, R2, R3, p 这些向量在公式里被当作了坐标值或坐标差来用，这非常 confusing。
        # 再次检查原始代码 sub_S = R1 + ...; sub_S = R2 - ... + p;
        # 如果 R1 = R1_value * A1 是从 yld_sit 出发的位移向量，那么点应该是 yld_sit + R1。
        # 同样，chj_sit + R2。
        # 这两个公式看起来不直接给出交点。它们像是在某个位移向量上加减一些项。
        # 但是 t_chj = norm(sub_S - chj_sit) 的计算又要求 sub_S 是一个点。
        # 鉴于此，我们**只能严格复制你的原始向量运算**来得到 sub_S，
        # 并假设这个结果 sub_S 就是计算出的点坐标。
        # 如果验证失败，很可能问题就出在这些公式的几何意义上。
        sub_S = R2 - (R2_value / R1_value) * (R2_value / (R1_value + R2_value)) * R3 + p

    # 最后一个检查，确保 sub_S 不是无穷大或 NaN (由于可能的除零等)
    if not np.all(np.isfinite(sub_S)):
        print("警告: 计算的 sub_S 包含非有限值 (inf 或 NaN)。")
        return None

    return sub_S


# ----------------------------------------------------------------------
# 实现标准的双站 AOA 定位 (用于对比)

def calculate_standard_aoa(yld_sit, chj_sit, yld_azimuth, yld_elevation, chj_azimuth, chj_elevation):
    """
    使用标准的双站 AOA 方法，计算两条射线最近点之间的中点作为估计位置。
    假设输入 azimuth 和 elevation 是标准方位角/仰角。
    """
    # 获取标准方向向量 (+X东, +Y北, +Z上)
    A1 = direction_vector_from_angles(yld_azimuth, yld_elevation)
    A2 = direction_vector_from_angles(chj_azimuth, chj_elevation)

    p = chj_sit - yld_sit

    # 标准射线方程: P1 = yld_sit + s1*A1, P2 = chj_sit + s2*A2
    # 我们想找到 s1 和 s2 使得 P2 - P1 垂直于 A1 和 A2
    # P2 - P1 = (chj_sit + s2*A2) - (yld_sit + s1*A1) = p + s2*A2 - s1*A1
    # (p + s2*A2 - s1*A1) . A1 = 0  => p.A1 + s2(A2.A1) - s1(A1.A1) = 0
    # (p + s2*A2 - s1*A1) . A2 = 0  => p.A2 + s2(A2.A2) - s1(A1.A2) = 0

    # A1, A2 是单位向量, A1.A1 = 1, A2.A2 = 1
    # 令 d = A1.A2
    # p.A1 + s2*d - s1 = 0   => s1 - d*s2 = p.A1
    # p.A2 + s2 - s1*d = 0   => -d*s1 + s2 = -p.A2

    A1_dot_p = np.dot(A1, p)
    A2_dot_p = np.dot(A2, p)
    A1_dot_A2 = np.dot(A1, A2)

    # 求解 2x2 线性系统:
    # [ 1  -d ] [s1] = [ p.A1 ]
    # [-d   1 ] [s2] = [-p.A2 ]

    matrix_s = np.array([[1, -A1_dot_A2], [-A1_dot_A2, 1]])
    vector_s = np.array([A1_dot_p, -A2_dot_p])

    det_s = np.linalg.det(matrix_s)

    if abs(det_s) < eps:
        # 矩阵接近奇异，射线可能平行或反平行
        print("警告: 标准方法: 射线接近平行，无法计算唯一解。")
        # 可以在此返回 None 或尝试其他方法 (如取任意一点)
        return None

    # 求解 s1 和 s2
    try:
        s1, s2 = np.linalg.solve(matrix_s, vector_s)
    except np.linalg.LinAlgError:
        print("警告: 标准方法: 线性系统求解失败。")
        return None

    # 最近点在射线上的位置
    P1_closest = yld_sit + s1 * A1
    P2_closest = chj_sit + s2 * A2

    # 估计位置取两个最近点的中点
    estimated_pos = (P1_closest + P2_closest) / 2.0

    # 也可以返回两个最近点之间的距离，作为射线不相交程度的指标
    # distance_between_rays = np.linalg.norm(P1_closest - P2_closest)

    return estimated_pos


# ----------------------------------------------------------------------
# 验证过程

def verify_localization_method(S_true):
    """
    对给定的真实闪电位置 S_true 运行验证过程。
    """
    print(f"\n--- 验证点: S_true = {S_true} ---")

    # 1. 模拟观测角度 (作为你的代码的输入)
    yld_az_obs, yld_el_obs = simulate_observation(S_true, yld_sit)
    chj_az_obs, chj_el_obs = simulate_observation(S_true, chj_sit)

    print(f"模拟观测角度 (yld): Az={yld_az_obs:.2f}°, El={yld_el_obs:.2f}°")
    print(f"模拟观测角度 (chj): Az={chj_az_obs:.2f}°, El={chj_el_obs:.2f}°")

    # 2. 使用你的代码的核心逻辑计算位置
    sub_S_user = calculate_sub_S_user_logic(yld_sit, chj_sit, yld_az_obs, yld_el_obs, chj_az_obs, chj_el_obs)

    # 3. 使用标准 AOA 方法计算位置 (用于对比)
    sub_S_standard = calculate_standard_aoa(yld_sit, chj_sit, yld_az_obs, yld_el_obs, chj_az_obs, chj_el_obs)

    # 4. 比较结果
    print("\n计算结果:")
    if sub_S_user is not None:
        print(f"你的方法计算位置: {sub_S_user}")
        error_user = np.linalg.norm(sub_S_user - S_true)
        print(f"你的方法误差 (与 S_true 的距离): {error_user:.2f} 米")
    else:
        print("你的方法计算失败。")

    if sub_S_standard is not None:
        print(f"标准方法计算位置: {sub_S_standard}")
        error_standard = np.linalg.norm(sub_S_standard - S_true)
        print(f"标准方法误差 (与 S_true 的距离): {error_standard:.2f} 米")
    else:
        print("标准方法计算失败。")

    # 比较两种方法之间的差异 (如果都计算成功)
    if sub_S_user is not None and sub_S_standard is not None:
        diff_methods = np.linalg.norm(sub_S_user - sub_S_standard)
        print(f"你的方法与标准方法结果差异: {diff_methods:.2f} 米")


# ----------------------------------------------------------------------
# 角度转换函数验证

def verify_angle_conversions():
    """
    验证 cart2sph_standard 和 direction_vector_from_angles 函数。
    """
    print("\n--- 角度转换函数验证 ---")

    # 测试向量 (+X东, +Y北, +Z上)
    # 1. 纯东方向
    v_east = np.array([1, 0, 0])
    az, el = cart2sph_standard(v_east)
    print(f"向量 [1, 0, 0] (东): Az={az:.2f}°, El={el:.2f}° (期望 Az=90°, El=0°)")  # arctan2(1,0) = 90deg
    v_reconstructed = direction_vector_from_angles(az, el)
    print(f"重构向量: {v_reconstructed} (期望接近 [1, 0, 0])")

    # 2. 纯北方向
    v_north = np.array([0, 1, 0])
    az, el = cart2sph_standard(v_north)
    print(f"向量 [0, 1, 0] (北): Az={az:.2f}°, El={el:.2f}° (期望 Az=0°, El=0°)")  # arctan2(0,1) = 0deg
    v_reconstructed = direction_vector_from_angles(az, el)
    print(f"重构向量: {v_reconstructed} (期望接近 [0, 1, 0])")

    # 3. 纯上方向
    v_up = np.array([0, 0, 1])
    # cart2sph_standard 会遇到 horizontal_distance = 0 的情况，arctan2(1, 0) 可以处理 z=1, horz=0
    # 但方位角 arctan2(x, y) = arctan2(0, 0) 是未定义的。
    # 约定正上方方位角为 0 或任意值，仰角为 90。
    az, el = cart2sph_standard(v_up)
    print(f"向量 [0, 0, 1] (上): Az={az:.2f}°, El={el:.2f}° (期望 Az=0°或任意值, El=90°)")
    # 重构时使用一个约定的方位角（如 0）
    v_reconstructed = direction_vector_from_angles(0, 90)
    print(f"重构向量 (Az=0, El=90): {v_reconstructed} (期望接近 [0, 0, 1])")

    # 4. 某个通用方向
    v_general = np.array([3, 4, 5])  # 3东, 4北, 5上
    az, el = cart2sph_standard(v_general)
    print(f"向量 [3, 4, 5]: Az={az:.2f}°, El={el:.2f}°")  # arctan2(3, 4) approx 36.87 deg
    v_reconstructed = direction_vector_from_angles(az, el)
    print(f"重构向量方向: {v_reconstructed}")
    # 比较重构向量方向和原始向量的单位向量方向
    print(f"原始单位向量方向: {v_general / np.linalg.norm(v_general)}")
    error_direction = np.linalg.norm(v_reconstructed - v_general / np.linalg.norm(v_general))
    print(f"重构方向误差: {error_direction:.2e}")


# ----------------------------------------------------------------------
# 运行验证

# 1. 验证角度转换函数
verify_angle_conversions()

# 2. 运行定位方法验证
# 选择一个真实的闪电位置 S_true (+X东, +Y北, +Z上)
S_true_test1 = np.array([5000, 2000, 8000], dtype=float)  # 距 yld_sit 较远的点
verify_localization_method(S_true_test1)

# 选择另一个不同的测试点
S_true_test2 = np.array([-3000, 10000, 6000], dtype=float)  # 另一个方向和距离的点
verify_localization_method(S_true_test2)

# 选择一个大致在 yld_sit 和 chj_sit 连线上的点 (可能更接近相交)
# chj_sit - yld_sit = [1991, -7841.2, 0]
# 选择一个点大约是 yld_sit + 0.5 * (chj_sit - yld_sit) + [0, 0, height]
mid_point_horizontal = yld_sit + 0.5 * (chj_sit - yld_sit)
S_true_test3 = mid_point_horizontal + np.array([0, 0, 7000], dtype=float)
verify_localization_method(S_true_test3)

# 选择一个可能导致射线接近平行的点 (例如，离站点非常远，且方向大致平行于站点连线)
# 如果站点在 XY 平面，选择一个很高的点，且其XY位置与站点连线方向相似
# 这可能难以构造，因为两个点确定一条线。要让观测方向平行，点需要在无穷远处。
# 考虑两个站点连线的方向向量 dir_stations = normalize(chj_sit - yld_sit)
# 选一个点 S_true 使得 (S_true - yld_sit) 和 (S_true - chj_sit) 方向向量相似
# 这意味着 S_true 离站点非常远，且大致在站点连线的延长线上。
dir_stations = (chj_sit - yld_sit) / np.linalg.norm(chj_sit - yld_sit)
S_true_test4 = yld_sit + dir_stations * 100000 + np.array([0, 0, 10000], dtype=float)  # 远距离点
verify_localization_method(S_true_test4)

# 考虑一个可能导致 R1_value <= R2_value 或 > 的点，并且 R1_value, R2_value 值差异大
# S_true_test1 的 R1_value 应该是 norm(S_true_test1 - yld_sit), R2_value 应该是 norm(S_true_test1 - chj_sit) (如果相交)
# 对于 test1: norm([5000, 2000, 8000]) approx 9.7k, norm([5000-1991, 2000-(-7841.2), 8000]) = norm([3009, 9841.2, 8000]) approx 13.2k
# R1_value < R2_value 在 test1 中应该发生
# 对于 test2: norm([-3000, 10000, 6000]) approx 12.0k, norm([-3000-1991, 10000-(-7841.2), 6000]) = norm([-4991, 17841.2, 6000]) approx 19.5k
# R1_value < R2_value 在 test2 中应该发生
# 构造 R1_value > R2_value 的点: S_true 更靠近 chj_sit
S_true_test5 = chj_sit + np.array([-1000, 500, 5000], dtype=float)  # 靠近 chj_sit 的点
# norm(S_true_test5 - yld_sit) = norm([1991-1000, -7841.2+500, 5000]) = norm([991, -7341.2, 5000]) approx 9k
# norm(S_true_test5 - chj_sit) = norm([-1000, 500, 5000]) approx 5.1k
# R1_value > R2_value 在 test5 中应该发生
verify_localization_method(S_true_test5)
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2025/10/4 23:19
# @Author : Ray
# @File : Analysis_Signal_Interval.py
# @Software: PyCharm
"""
分析正负先导的发展特征
"""
import pandas as pd
import numpy as np

# 1. 加载2d定位结果数据
result = pd.read_csv('result_yld_3.8e8_4e8_window_512_128_去零飘_滤波_加窗_阈值15_30_80.txt', sep=r'\s+')

# 2. 对2d结果进行过滤筛选
filtered_result = result[
    (abs(result['t123']) < 1) &  # t123的绝对值小于1
    (abs(result['Rcorr']) > 0.65) &  # Rcorr的绝对值大于0.65
    (result['Start_loc'] < 4e8) &  # Start_loc小于4e8
    (result['Start_loc'] > 3.8e8)  # Start_loc大于3.8e8
]

# 3. 分析负先导的发展特征
# 选出负先导数据
negative_lead_result = filtered_result[
    (filtered_result.Start_loc >= 3.9e8) &
    (filtered_result.Start_loc < 4e8)
]
time_diffs = np.diff(negative_lead_result.Start_loc)
negative_time_diffs = pd.DataFrame(data=time_diffs)
# 查看数据的分布
print("负先导阶段的每个结果之间的采样点之差:\n", negative_time_diffs.describe())


# 4. 分析正先导的发展特征
# 选出正先导数据
positive_lead_result = filtered_result[
    (filtered_result.Start_loc >= 3.8e8) &
    (filtered_result.Start_loc < 3.9e8)
]
time_diffs = np.diff(positive_lead_result.Start_loc)
positive_time_diffs = pd.DataFrame(data=time_diffs)
# 查看数据的分布
print("正先导阶段的每个结果之间的采样点之差:\n", positive_time_diffs.describe())


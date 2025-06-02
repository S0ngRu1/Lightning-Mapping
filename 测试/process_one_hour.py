#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
@Project ：Data-Factory 
@File    ：process_one_hour.py
@IDE     ：PyCharm 
@Author  ：Cai Songrui 
@Date    ：2025/3/3 21:09 
'''

import pandas as pd
from matplotlib import pyplot as plt


def process_one_hour(file_path):
    df = pd.read_excel(file_path).to_dict("list")
    time_list = df["第一列"]
    result = []
    filtered_time_list = []
    for i in range(len(time_list)):
        time = float(time_list[i].split(".")[1])
        if abs(time) < 1000:
            result.append(time)
            filtered_time_list.append(time_list[i])
        if abs(time - 1e9) < 1000:
            result.append(time - 1e9)
            filtered_time_list.append(time_list[i])
    # 保存过滤后的时间差 和时间戳在一个excel中
    df = pd.DataFrame({"时间戳": filtered_time_list, "时间差": result})
    df.to_excel("one_hour_filtered.xlsx", index=False)
    print("Save successfully")

if __name__ == '__main__':
    file_path = "one_hour.xlsx"
    process_one_hour(file_path)
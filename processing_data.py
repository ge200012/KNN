# import pandas as pd
# import os

# import pandas as pd

# #####################################################################################
# import signal
# import sys 
    
# #signal.siginterrupt(signal.SIGTERM,False)

# def sigterm_handler(signum,frame):
#     pass
#     #print("Received SIGTERM signal. Exiting gracefully.")

# signal.signal(signal.SIGTERM,sigterm_handler)
# #####################################################################################
   
# # 读取纯文本文件
# file_path = '/home/KNN/data/cp_data/chengdu/trajectory_data'  # 替换为你的文件路径
# trajectories = []

# with open(file_path, 'r') as file:
#     lines = file.readlines()

# # 处理每条轨迹
# i = 0
# while i < len(lines):
#     index_id = lines[i].strip()  # 读取索引ID
#     i += 1
    
#     trajectory = []
    
#     while i < len(lines) and lines[i].strip():  # 读取轨迹点
#         parts = lines[i].strip().split()
#         if len(parts) == 3:  # 确保有三列
#             latitude = float(parts[1])
#             longitude = float(parts[2])
#             # 添加轨迹点
#             trajectory.append([latitude, longitude])
#         i += 1
    
#     # 去重连续相同的轨迹点
#     cleaned_trajectory = []
#     for point in trajectory:
#         if not cleaned_trajectory or cleaned_trajectory[-1] != point:
#             cleaned_trajectory.append(point)

#     # 过滤轨迹长度在90到300之间的轨迹
#     if 90 <= len(cleaned_trajectory) <= 300:
#         trajectories.append(cleaned_trajectory)

# # 创建 DataFrame
# print(len(trajectories))
# for i in range(min(30,len(trajectories))):
#     print(trajectories[i])
# t_table = pd.DataFrame({'tra': trajectories})

# # 保存为 pickle 文件
# t_table.to_pickle('/home/KNN/data/cp_data/chengdu/scale_t/processed_trajectories.pkl')

# # 打印前 3000 行
# print(t_table.head(3000))


# import pandas as pd
# import os
# # def process_trajectories(file_path, output_path):
# #     trajectories = []
# #     current_trajectory = []
    
# #     with open(file_path, 'r') as f:
# #         for line in f:
# #             line = line.strip()
# #             if line.isdigit():  # 轨迹 ID 行
# #                 if current_trajectory:  # 如果当前轨迹非空，处理并存储
# #                     # 去重
# #                     unique_trajectory = []
# #                     for point in current_trajectory:
# #                         if not unique_trajectory or unique_trajectory[-1] != point:
# #                             unique_trajectory.append(point)
# #                     # 只保留长度在 90 到 300 之间的轨迹
# #                     if 90 <= len(unique_trajectory) <= 300:
# #                         trajectories.append(unique_trajectory)
# #                     current_trajectory = []  # 重置当前轨迹
# #             else:  # 轨迹点行
# #                 parts = line.split()
# #                 if len(parts) == 3:
# #                     time, lat, lon = parts
# #                     current_trajectory.append([float(lat), float(lon)])  # 只保留经纬度

# #     # 处理最后一个轨迹
# #     if current_trajectory:
# #         unique_trajectory = []
# #         for point in current_trajectory:
# #             if not unique_trajectory or unique_trajectory[-1] != point:
# #                 unique_trajectory.append(point)
# #         if 90 <= len(unique_trajectory) <= 300:
# #             trajectories.append(unique_trajectory)

# #     # 保存为 DataFrame 并 pickle
# #     t_table = pd.DataFrame({'tra': trajectories})
# #     t_table.to_pickle(output_path)
    
# # #####################################################################################
# import signal
# import sys 
    
# # #signal.siginterrupt(signal.SIGTERM,False)

# def sigterm_handler(signum,frame):
#     pass
# #     #print("Received SIGTERM signal. Exiting gracefully.")

# signal.signal(signal.SIGTERM,sigterm_handler)
# #####################################################################################
# # # 读取纯文本文件
# # file_path = '/home/KNN/data/cp_data/chengdu/30_90_lcss'  # 替换为你的文件路径
# # trajectories = []

# # with open(file_path, 'r') as file:
# #     lines = file.readlines()


# # for i in range(30):
# #     print(lines[i])
# # exit(1)

# # 使用示例
# #input_file_path = '/home/KNN/data/cp_data/xian/scale_t/trajectory_data'
# output_file_path = '/home/KNN/data/cp_data/chengdu/scale_t/processed_trajectories.pkl'
# #process_trajectories(input_file_path, output_file_path)

# #读取并打印处理后的数据
# t_table = pd.read_pickle(output_file_path)
# print(len(t_table))
# for i in range(30):
#     inner_list=t_table.iloc[i]
#     for j in range(len(inner_list)):
#         for t in range(len(inner_list[j])):
#             print(inner_list[j][t])
# print(t_table.head(3))


import pandas as pd
# #####################################################################################
import signal
import sys 
    
# #signal.siginterrupt(signal.SIGTERM,False)

def sigterm_handler(signum,frame):
    pass
#     #print("Received SIGTERM signal. Exiting gracefully.")

signal.signal(signal.SIGTERM,sigterm_handler)
#####################################################################################
# 从处理后的数据文件中读取数据
input_file_path = '/home/KNN/data/cp_data/xian/scale_t/processed_trajectories.pkl'  # 输入文件路径
output_file_path = '/home/KNN/data/cp_data/xian/scale_t/processed_trajectories_N.pkl'  # 输出文件路径

# 读取处理后的数据
t_table = pd.read_pickle(input_file_path)

# 归一化经纬度的函数
def normalize_trajectories(df):
    all_latitudes = []
    all_longitudes = []

    # 收集所有经纬度
    for trajectory in df['tra']:
        for point in trajectory:
            all_latitudes.append(point[0])
            all_longitudes.append(point[1])

    # 计算经纬度的最小值和最大值
    min_latitude = min(all_latitudes)
    max_latitude = max(all_latitudes)
    min_longitude = min(all_longitudes)
    max_longitude = max(all_longitudes)

    # 归一化
    normalized_trajectories = []
    for trajectory in df['tra']:
        normalized_trajectory = []
        for point in trajectory:
            normalized_latitude = (point[0] - min_latitude) / (max_latitude - min_latitude)
            normalized_longitude = (point[1] - min_longitude) / (max_longitude - min_longitude)
            normalized_trajectory.append([normalized_latitude, normalized_longitude])
        normalized_trajectories.append(normalized_trajectory)

    return pd.DataFrame({'tra': normalized_trajectories})

# 对轨迹进行归一化
normalized_table = normalize_trajectories(t_table)

# 保存归一化后的数据到新的 pickle 文件
normalized_table.to_pickle(output_file_path)

# 打印归一化后的数据的前 10 行
print(normalized_table.head(10))

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# 计算 cross ratios
def calculate_cross_ratios(data, settings):
    ratios = {}
    for graph_id, metrics in data.items():
        baseline_index = metrics['settings'].index('1')
        baseline_internal_cross = metrics['internal cross'][baseline_index]
        baseline_external_cross = metrics['external cross'][baseline_index]

        indexs = [metrics['settings'].index(setting) for setting in settings]

        if baseline_internal_cross == 0 or baseline_external_cross == 0:
            continue  # 跳过任何 baseline cross 为0的情况

        metrics['internal cross'] = np.array(metrics['internal cross'])[indexs]
        metrics['external cross'] = np.array(metrics['external cross'])[indexs]

        graph_ratios = {
            'internal cross ratio': np.array(metrics['internal cross']) / baseline_internal_cross,
            'external cross ratio': np.array(metrics['external cross']) / baseline_external_cross
        }
        ratios[graph_id] = graph_ratios
    return ratios

# 统计全局 cross ratios 平均值
def global_cross_ratio_average(ratios):
    all_internal_ratios = []
    all_external_ratios = []
    for graph_id, graph_ratios in ratios.items():
        all_internal_ratios.append(graph_ratios['internal cross ratio'])
        all_external_ratios.append(graph_ratios['external cross ratio'])
    
    mean_internal_ratio = np.mean(all_internal_ratios, axis=0)
    mean_external_ratio = np.mean(all_external_ratios, axis=0)
    return mean_internal_ratio, mean_external_ratio

# 主函数
def main():
    with open('all_metrics.json', 'r') as file:
        data = json.load(file)
    print('data length:', len(data))

    original_settings = data[next(iter(data))]['settings']  # 假设每个图的设置相同
    settings = [setting for setting in original_settings if not setting.startswith('t')]

    print('Settings:', original_settings, settings)
    ratios = calculate_cross_ratios(data, settings)

    # print('Cross Ratios:', ratios)
    mean_internal_ratio, mean_external_ratio = global_cross_ratio_average(ratios)

    print('Internal Ratios:', mean_internal_ratio)
    print('External Ratios:', mean_external_ratio)
    # 绘制全局曲线图
    plt.figure(figsize=(10, 5))
    plt.plot(settings, mean_internal_ratio, label='Average Internal Cross Ratio')
    plt.plot(settings, mean_external_ratio, label='Average External Cross Ratio')
    plt.xlabel('Focus Edge Weight')
    plt.ylabel('Cross Ratio')
    plt.title('Global Cross Ratios by Focus Edge Weight')
    plt.legend()
    plt.grid(True)
    plt.xticks(settings)  # 设置 x 轴刻度
    plt.show()
    plt.savefig('weight.png')

if __name__ == '__main__':
    main()

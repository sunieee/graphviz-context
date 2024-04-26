import matplotlib.pyplot as plt
import os
import re
import xml.etree.ElementTree as ET


def parse_filename(filename):
    """ 解析文件名以便进行正确的排序。
        首先根据前面的主数字排序，如果存在，则根据后面的次数字排序。
    """
    match = re.match(r'(\d+)(?:_(\d+))?\.log', filename)
    if match:
        main_num = int(match.group(1))
        sub_num = int(match.group(2)) if match.group(2) else 0  # 如果没有第二个数字，默认为0
        return (main_num, sub_num)
    return (float('inf'),)  # 如果不符合数字_数字.log格式，放在排序最后

def extract_metrics(filename):
    """ 从文件中提取指标 """
    metrics = {}
    with open(filename, 'r') as file:
        for line in file:
            if 'all cross:' in line:
                metrics['all cross'] = int(line.split(':')[-1].strip())
            elif 'weigted cross:' in line:
                metrics['weigted cross'] = int(line.split(':')[-1].strip())
            elif 'internal cross:' in line:
                metrics['internal cross'] = int(line.split(':')[-1].strip())
            elif 'external cross:' in line:
                metrics['external cross'] = int(line.split(':')[-1].strip())
            elif 'time cost:' in line:
                metrics['time cost'] = float(line.split(':')[-1].strip())
    return metrics

def draw(_dir):
    # 列出当前目录下所有.log文件并排序
    files = [f for f in os.listdir(_dir) if f.endswith('.log')]
    files.sort(key=parse_filename)
    svgs = [f.replace('.log', '.svg') for f in files]

    # 提取所有文件的指标
    data = {'all cross': [], 'weigted cross': [], 'internal cross': [], 'external cross': [],
            'time cost': []}
    filenames = []
    
    for file in files:
        file = os.path.join(_dir, file)
        metrics = extract_metrics(file)
        for key in data:
            data[key].append(metrics.get(key, 0))
        filenames.append(file.split('/')[-1].split('.')[0].replace('1_', '('))

    data['width'] = []
    data['height'] = []

    for svg in svgs:
        # 获取svg文件的宽度和高度
        tree = ET.parse(os.path.join(_dir, svg))
        root = tree.getroot()
        width = int(root.attrib['width'].replace('px', '').replace('pt', ''))
        height = int(root.attrib['height'].replace('px', '').replace('pt', ''))
        data['width'].append(width)
        data['height'].append(height)


    # 为每个指标绘图并保存
    color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    def plot(*keys):
        plt.figure()
        for ix, key in enumerate(keys):
            plt.plot(filenames, data[key], color=color_list[ix], marker='o', label=key)
        plt.xlabel('Filename')
        plt.ylabel(key)
        # plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{_dir}/{keys[0]}.png")
        plt.close()

    plot('all cross')
    plot('weigted cross')
    plot('internal cross')
    plot('external cross')
    plot('width', 'height')
    plot('time cost')


if __name__ == '__main__':
    # 枚举当前目录下的所有文件夹
    dirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    for _dir in dirs:
        print(f"Processing {_dir}...")
        try:
            draw(_dir)
        except Exception as e:
            print(f"Error: {e}")

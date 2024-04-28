#!/bin/bash

# http://ye-sun.com:1401/sy/2024/graphviz-9.0.0/example/

# 获取脚本所在的绝对路径
script_dir=$(dirname "$(readlink -f "$0")")

# 切换到脚本所在目录
cd "$script_dir"

for d in */ ; do 
    cd "$d"
    echo "Working in directory: $(pwd)"
    

    echo "  Creating .dot files..."
    prototype="0.dot"
    # 检查原型文件是否存在
    if [ ! -f "$prototype" ]; then
        echo "Error: Prototype file $prototype does not exist."
        exit 1
    fi
    # 删除0.dot以外的所有文件
    # rm -f !("0.dot")
    # ./run.sh: line 22: syntax error near unexpected token `('
    for file in *; do
        if [ "$file" != "0.dot" ]; then
            rm -f "$file"
        fi
    done


    # 生成带有不同权重的文件
    for weight in 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 35 40; do
        sed -e 's/crossing_type=0/crossing_type=1/g' -e "s/edge \[weight=10\]/edge \[weight=$weight\]/g" $prototype > "1_${weight}.dot"
    done
    # 生成其他crossing_type的版本
    for type in 2 3; do
        sed "s/crossing_type=0/crossing_type=$type/g" $prototype > "${type}.dot"
    done

    # 遍历当前目录下所有 .dot 文件
    for file in *.dot; do
        # 检查是否真的存在 .dot 文件，避免无文件时的错误执行
        if [[ -f "$file" ]]; then
            # 获取不带扩展名的文件名
            basename="${file%.*}"
            echo "  processing ${d}${file}"
            # 获取系统秒数，精确到毫秒
            timestamp=$(date +%s%3N)
            # 使用 dot 命令将 .dot 文件转换为 .svg，并将错误信息记录到 .log 文件
            cat "$file" | dot -v -Tsvg > "${basename}.svg" 2> "${basename}.log"
            # 追加毫秒数到 .log 文件
            cost=$(( $(date +%s%3N) - $timestamp ))
            echo "time cost: ${cost}" >> "${basename}.log"
        else
            echo "No .dot files found in the script's directory."
            break
        fi
    done
    cd ..
done

python draw_metrics.py


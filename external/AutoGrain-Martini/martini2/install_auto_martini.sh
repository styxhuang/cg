#!/bin/bash

# 创建虚拟环境
python -m venv auto_martini_env

# 激活虚拟环境
source auto_martini_env/bin/activate

# 安装必要的依赖
pip install numpy
pip install MDAnalysis
pip install networkx
pip install tqdm
pip install matplotlib

# 克隆 auto_martini 仓库
git clone https://github.com/ajz34/auto_martini.git

# 进入 auto_martini 目录
cd auto_martini

# 安装 auto_martini
pip install -e .

echo "安装完成！" 
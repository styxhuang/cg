# 项目描述

我现在想实现一个自动粗粒化的脚本

输入：
    --aa: pdb 文件
    --type: 纯蛋白/改性残基/蛋白+配体
    --smi: SMILES结构  
输出：
    --cg: 粗粒化之后的pdb文件
    --top: 粗粒化之后的力场文件
其他：
    --verbose, -v: 更详细的输出

流程：
1. 如果是纯蛋白结构，使用martinize.py生成粗粒化力场
2. 如果是改性残基，或者是SMILES，则切片，然后各个切片使用auto_martini进行粗粒化，最终把生成的结构，力场合并起来
3. 如果是蛋白+配体：
    1. 默认配体的残基名称为LIG，切片，然后用auto_martini粗粒化
    2. 纯蛋白部分，用matinize.py粗粒化


# 测试martinize.py 升级到python3支持

`python martinize/martinize.py -f assets/pure.pdb -o tests/cg.gro`

`martinize v2: python martinize/martinize.py -f test/ARG/ARG.pdb -o tmp/tmp.top`
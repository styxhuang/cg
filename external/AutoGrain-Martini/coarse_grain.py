#!/usr/bin/env python3
import argparse
import os
import subprocess
from typing import Optional, Tuple

def parse_args():
    parser = argparse.ArgumentParser(description='自动粗粒化脚本')
    parser.add_argument('--aa', required=True, help='输入的PDB文件路径')
    parser.add_argument('--type', required=True, choices=['protein', 'modified', 'protein_ligand'],
                      help='分子类型：纯蛋白/改性残基/蛋白+配体')
    parser.add_argument('--smi', help='SMILES结构（如果需要）')
    parser.add_argument('--cg', required=True, help='输出的粗粒化PDB文件路径')
    parser.add_argument('--top', required=True, help='输出的粗粒化力场文件路径')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示详细输出')
    return parser.parse_args()

def run_martinize(input_pdb: str, output_pdb: str, output_top: str, verbose: bool = False) -> None:
    """运行martinize.py进行粗粒化
    
    参数说明：
    - input_pdb: 输入的PDB文件
    - output_pdb: 输出的粗粒化PDB文件
    - output_top: 输出的拓扑文件
    - verbose: 是否显示详细输出
    """
    # martinize.py的基本命令
    cmd = [
        'python', 'martinize/martinize.py',
        '-f', input_pdb,  # 输入PDB文件
        '-o', output_pdb,  # 输出PDB文件
        '-x', output_top,  # 输出拓扑文件
        '-p', 'topol.top',  # 生成GROMACS拓扑文件
        '-cys', 'auto',  # 自动处理二硫键
        '-elastic',  # 添加弹性网络
        '-ef', '500',  # 弹性网络力常数
        '-el', '0.5',  # 弹性网络截断距离
        '-eu', '0.9',  # 弹性网络上界
        '-ea', '0',  # 弹性网络角度
        '-ep', '0'  # 弹性网络平面
    ]
    
    if verbose:
        print(f"执行命令: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"martinize.py执行失败: {e}")
        raise

def run_auto_martini(input_pdb: str, output_pdb: str, output_top: str, verbose: bool = False) -> None:
    """运行auto_martini进行粗粒化"""
    cmd = ['auto_martini', '-f', input_pdb, '-o', output_pdb, '-x', output_top]
    if verbose:
        print(f"执行命令: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def process_protein(input_pdb: str, output_pdb: str, output_top: str, verbose: bool = False) -> None:
    """处理纯蛋白结构"""
    if verbose:
        print("处理纯蛋白结构...")
    run_martinize(input_pdb, output_pdb, output_top, verbose)

def process_modified(input_pdb: str, output_pdb: str, output_top: str, verbose: bool = False) -> None:
    """处理改性残基结构"""
    if verbose:
        print("处理改性残基结构...")
    run_auto_martini(input_pdb, output_pdb, output_top, verbose)

def process_protein_ligand(input_pdb: str, output_pdb: str, output_top: str, verbose: bool = False) -> None:
    """处理蛋白+配体结构"""
    if verbose:
        print("处理蛋白+配体结构...")
    
    # 分离蛋白和配体
    protein_pdb = "protein.pdb"
    ligand_pdb = "ligand.pdb"
    
    # 使用grep分离LIG残基
    subprocess.run(f"grep -v 'LIG' {input_pdb} > {protein_pdb}", shell=True, check=True)
    subprocess.run(f"grep 'LIG' {input_pdb} > {ligand_pdb}", shell=True, check=True)
    
    # 分别处理蛋白和配体
    protein_cg = "protein_cg.pdb"
    protein_top = "protein_cg.top"
    ligand_cg = "ligand_cg.pdb"
    ligand_top = "ligand_cg.top"
    
    run_martinize(protein_pdb, protein_cg, protein_top, verbose)
    run_auto_martini(ligand_pdb, ligand_cg, ligand_top, verbose)
    
    # 合并结果
    subprocess.run(f"cat {protein_cg} {ligand_cg} > {output_pdb}", shell=True, check=True)
    subprocess.run(f"cat {protein_top} {ligand_top} > {output_top}", shell=True, check=True)
    
    # 清理临时文件
    for f in [protein_pdb, ligand_pdb, protein_cg, protein_top, ligand_cg, ligand_top]:
        os.remove(f)

def main():
    args = parse_args()
    
    if args.verbose:
        print(f"输入文件: {args.aa}")
        print(f"分子类型: {args.type}")
        if args.smi:
            print(f"SMILES: {args.smi}")
        print(f"输出PDB: {args.cg}")
        print(f"输出拓扑: {args.top}")
    
    if args.type == 'protein':
        process_protein(args.aa, args.cg, args.top, args.verbose)
    elif args.type == 'modified':
        process_modified(args.aa, args.cg, args.top, args.verbose)
    elif args.type == 'protein_ligand':
        process_protein_ligand(args.aa, args.cg, args.top, args.verbose)

if __name__ == '__main__':
    main() 
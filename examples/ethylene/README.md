# 乙烯粗粒化示例

## 文件
- `mapping_two_beads.xml`：两珠映射，将乙烯分为两个 CH2 bead（A、B）。

## 前提
- 已安装 VOTCA（Homebrew 安装或 Conda 环境 `cg`），命令 `csg_map` 可用。
- 准备好原子级拓扑与轨迹：`topol.tpr`（或 `*.top/*.gro`）与 `traj.xtc`（或 `*.trr`）。
- 原子命名需与映射文件一致：`C1 C2 H1 H2 H3 H4` 属于同一分子。

## 运行
### 直接命令
```bash
csg_map --top topol.tpr --trj traj.xtc --cg example/ethylene/mapping_two_beads.xml
```

### Conda 环境
```bash
conda run -n cg bash -lc "csg_map --top topol.tpr --trj traj.xtc --cg example/ethylene/mapping_two_beads.xml"
```

### 后端接口
1. 在页面上传 `topol.tpr / traj.xtc / mapping_two_beads.xml`。
2. 调用 `POST /api/votca/map`，传递上传返回的 `rid`。

## 提示
- 可用 `csg_dump --top topol.tpr --cg example/ethylene/mapping_two_beads.xml` 检查映射与命名。
- 若命名不同，请修改 `mapping_two_beads.xml` 中 `<atom name="..."/>` 为你的拓扑命名。
- 若需单珠映射，可将所有原子放入同一 `<bead>`。
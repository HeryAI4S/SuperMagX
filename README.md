# SuperMagX

**语言 / Languages：** **简体中文（当前）** | [English](README.en.md)

SuperMagX 是一个面向**脉冲磁体计算**的模块化 Python 项目。第一阶段已实现参数化的二维轴对称**静磁场**求解器，以及独立的输运电流损耗与感应损耗模型。热分析和结构力学目前仅预留模块边界，**尚未实现温度场、位移场或应力场求解器**。

## 安装与运行

```powershell
python -m pip install -e ".[test]"
python -m pytest -q
python -m supermagx examples/two_zone_case.json result.json
```

程序支持任意数量、带名称的独立供电分区。每匝实体导体具有独立的矩形金属截面、轴向位置、供电分区及电流绕向。

- `radial_stack()` 按可调节的径向层数、每层匝数、尺寸、节距和分区间隙生成规则绕组。
- 对于不规则绕组，可使用 `Magnet(tuple(RectTurn(...), ...))`，或在 JSON 中指定 `geometry.kind = "explicit_turns"`。
- 相同的分区标识符代表多匝共用串联电流；不同分区可使用不同的分段线性电流波形。

## 模块与数据流

```text
几何 geometry + 激励 excitation
       |               |
       v               v
电磁 electromagnetics ---> 损耗 losses
       |                       |
       v                       v
 力学 mechanics            热 thermal
   （待实现）               （待实现）
```

| 模块 | 第一阶段的功能 |
| --- | --- |
| `geometry` | 有限矩形导体匝、任意径向组及供电分区映射 |
| `excitation` | 分段线性时间网格上的分区独立电流 |
| `electromagnetics` | 自由空间轴对称磁场（导体外观察点）及有限截面导体分区互感矩阵 |
| `losses.joule` | 固定电阻率下 `I²R` 的精确分段时间积分；电阻与选定截面电流密度分布一致 |
| `losses.eddy` | 导体匝中心磁场筛选：经标定的单极点响应，或无拟合参数的理想平板扩散模态；均含脉冲结束后的感应电流衰减热 |
| `losses.network` | 耦合导电单元磁扩散与脉冲发热，包含涡流尾部；可预编译算子用于重复波形/电阻率扫描 |
| `thermal`、`mechanics` | 保留接口，尚不能计算局部温度、位移或真实峰值应力 |

命令 `python -m supermagx` 读取 JSON 算例并生成 JSON 结果，其中包含磁场探针所用电流、积分阶数与阶间变化、分区损耗、未求解模块的状态。

- [双分区示例](examples/two_zone_case.json)
- [耦合损耗示例](examples/coupled_loss_case.json)

## 数值模型与精度

### 电流密度模型

每匝可在 `RectTurn` 或 `RadialGroup` 中指定 `source_profile`：

| 模式 | 分布 | 物理解释 |
| --- | --- | --- |
| `"dc_equilibrium"`（**代码默认**） | `J_phi = I/(r H)`，`H = ∫_metal dr dz / r` | 均匀电阻率闭合环导体的直流电阻平衡输运电流 |
| `"uniform"` | `J_phi = I/A` | 金属截面上均匀周向电流密度，例如对应 FEM 中给定均匀 `External Current Density` 的静磁算例 |

这里的 `A` 是金属 `r–z` 截面积，`H` 为截面上的 `1/r` 积分。源分布选择同时影响**磁场源项与焦耳等效电阻**。

**重要区别：**两种模式都是规定的静磁/输运电流分布；在真实脉冲中，还会产生磁扩散和邻近效应导致的感应电流。不能把其中任意一种静态分布直接等同于全时间的真实总电流密度。

### 磁场与互感积分

电磁求解器把轴对称圆形细丝 Green 核积分到每匝有限矩形金属截面上：

- `field()` 逐次加密 Gauss 积分阶数，直到磁场变化满足设置的容差，否则抛出 `ConvergenceError`。
- 互感采用有限截面导体对的自适应积分；自感的可积对数奇异性使用 **Duffy 变换**处理。
- 利用轴向坐标差的梯形权重，对两个轴向积分作解析降维。
- 近重合圆环使用由真实几何距离计算的椭圆积分互补参数，降低 `1-m` 相减的精度损失。
- 轴线附近使用 Maxwell 级数极限计算径向磁场。
- 程序**没有使用经验磁场放大系数**。

### 感应损耗：理想平板模态及单极点近似

以下 JSON 选择不使用拟合参数的理想平板扩散模态展开：

```json
"eddy": {"model": "slab_modes", "modes": 16}
```

使用原有 `"eddy": {"gamma": ...}` 则选择经过标定的单极点筛选模型。

平板模型分别处理导体两个截面方向的奇数扩散模态，但依然使用**导体中心位置的外加磁场**，并假设两个平板方向独立。波形最后一个节点之后的感应热衰减尾部也被积分计入。详见[第一性原理推导](docs/first_principles.md)。

### 耦合导电单元磁扩散

下列设置选择耦合单元模型：

```json
"eddy": {
  "model": "coupled_network",
  "radial_cells": 4,
  "axial_cells": 4,
  "cell_distribution": "edge_cosine",
  "refinement_meshes": [[2, 2], [3, 3], [4, 4]]
}
```

输运电流与感应电流使用相同的 `resistivity_ohm_m`。结果区分 `joule_j`（输运焦耳热）、`eddy_before_tail_j`（波形结束前的涡流热）、`eddy_tail_j`（结束后衰减尾部热）和 `total_j`（总热量）；`network_refinement` 报告网格间的变化。

由于感应电流模态可能跨越不同供电分区耦合，程序**不将全部感应热强行归属到单个供电分区**。应加密导电单元直到**感应热自身**收敛，而不应只查看总热量。

耦合网络需要构造稠密互感矩阵和求解特征值问题，大规模参数扫描可能较昂贵；平板近似仍适合初步筛选。同一几何下可先编译一次算子，再重复计算不同波形和电阻率：

```python
operator = compile_coupled_loss(magnet, ...)
result = operator.evaluate(waveform, resistivity=rho)
```

参见[预编译扫描示例](examples/compiled_loss_scan.py)。`edge_cosine` 将更多单元布置在金属边界附近，仍必须检查网格收敛。

**接口限制：**耦合网络要求每匝 `source_profile = "dc_equilibrium"`；若指定 `uniform`，程序会拒绝运行。这与其直流电阻平衡基底和输运/感应热精确分解有关。

连续理论推导、单环与双环 COMSOL 6.4 对照见[损耗理论与验证报告](docs/loss_theory_validation.md)。

## 适用范围与验证状态

这些数值收敛结果仅适用于**当前声明的自由空间、给定直流输运电流分布模型**。积分阶数变化满足容差并非严格误差上界，也不能保证所有算例与 COMSOL 的偏差小于 1%。

- 精度受控的 `field()` 仅允许**导体外**的观察点；涡流筛选使用的 `screening_coefficients()` 在匝中心取样，不具有相同的精度保证。
- 当前阶段未建模铁磁材料、引线、接头、真实电源回路、导电支撑、温度相关电阻率及机械预紧。
- 单极点与平板感应损耗模型属于筛选近似。
- 耦合单元模型已有单环/双环 COMSOL 6.4 独立检查，但不能据此认定完整 APEX 绕组和所有设计变型均已验证。
- 当前命令行磁场探针仍使用规定源分布的**静磁**结果；耦合磁扩散损耗模型尚未将逐单元瞬态电流完整回传到探针磁场及洛伦兹体积力。

工程使用前请阅读[精度与验证说明](docs/accuracy.md)。

## 当前 APEX 物理模型与约束

[APEX 物理模型规范](docs/apex/README.md) 登记了 **20 T 中心场 / 200 mm 室温净孔**的 Cu–Nb 脉冲磁体设计要求（包括**每路电源峰值电流严格小于 25 kA**）、M1 控制方程与代码映射，以及后续磁场、电感、磁能和电磁力的有限元验证门槛。

M1 复用现有 Green 核求解器。若需与 COMSOL 给定均匀 `External Current Density` 的静磁场模型进行同假设比较，应选择 `source_profile = "uniform"`；当前磁扩散损耗网络则要求 `dc_equilibrium`。两个模块不应在源分布不一致的情况下直接混合为已验证的全瞬态解。

## APEX Cu–Nb 窄截面历史研究（OP5）

[OP5 筛选研究](studies/apex_cunb_op5/op5_report.md) 考察 70 种径向金属宽度至少为 5 mm 的 Cu–Nb 矩形导体方案，优化对称脉冲斜坡时间，比较单次总热量、电磁力代理指标、分区电压和逐匝平均温升。

[可复现扫描脚本](studies/apex_cunb_op5/scan_narrow.py) 采用 OP3 标定的单极点损耗近似，结果表格保存在相同目录。该研究未对新几何完成结构或瞬态热 FEM 验证；其历史推荐点使用 **29.785 kA**，**不满足目前严格小于 25 kA 的电流约束**，不得误认为当前可行方案。

## 许可证

MIT 许可证；详见 [LICENSE](LICENSE)。

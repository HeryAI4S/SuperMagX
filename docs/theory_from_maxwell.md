# 从 Maxwell 方程到 SuperMagX 解析与半解析模型

本文固定电流、几何、材料、损耗及验证的定义，供 APEX CuNb 脉冲磁体的参数筛选使用。所有公式采用 SI 单位。模型中的“匝”是绝缘隔开的完整环向导体，截面为 (r-z) 平面中的金属区域 (A_n)，回转后才是三维导体。多个匝可在一个电源区内串联；不同区的电流波形可以独立指定。当前 Python 包已经实现参数化几何、自由空间磁静态场与互感、精确固定电阻率的传输电流焦耳热，以及三档涡流筛选模型。空间温度场与真实应力场尚未在 Python 包中实现。

| 模块 | 当前代码入口 | 已实现的内容与边界 |
| --- | --- | --- |
| 几何与电源 | `src/supermagx/geometry.py`、`src/supermagx/excitation.py` | 任意数量矩形匝、电源区和分段线性区电流；真实引线及开关电路未解 |
| 电磁 | `src/supermagx/electromagnetics/` | 自由空间环核、有限截面磁场与互感；时变电流重排交给损耗网络 |
| 损耗 | `src/supermagx/losses/` | 精确固定电阻率的传输热、单极点/理想薄板筛选、耦合胞元磁扩散及关流后尾热 |
| 传热与结构 | `src/supermagx/thermal/`、`src/supermagx/mechanics/` | 目前为模块接口；空间温度场和真实应力须由外部 FEM 求解并核对 |

## 1. Maxwell 方程、边界与适用范围

从微分形式出发，

\[
\nabla\cdot\boldsymbol B=0,\quad
\nabla\times\boldsymbol E=-\partial_t\boldsymbol B,\quad
\nabla\times\boldsymbol H=\boldsymbol J+\partial_t\boldsymbol D,\quad
\nabla\cdot\boldsymbol D=\rho_e .
\]

本阶段采用磁准静态近似：脉冲尺度相对于电磁传播时间足够慢，忽略位移电流；导体各向同性且 (oldsymbol B=\mu_0\boldsymbol H)、(oldsymbol E=\rho\boldsymbol J)。这里的 \(\rho\) 是**电阻率**，与上式电荷密度 \(\rho_e\) 不同。CuNb 基准值为 \(6\times10^{-9}\ \Omega\,\mathrm m\)（77 K）；不把该固定值解释成已验证的 \(\rho(T,B)\) 关系。自由空间、无铁磁体、无导电支撑件和理想绝缘匝，是下列简化的物理边界。若引线、横向导电连接、金属约束层或材料非线性显著，需要相应扩展场和电路模型。

轴对称且仅有环向传输与感应电流时，令 \(\boldsymbol A=A_\phi(r,z,t)\boldsymbol e_\phi\)，则

\[
B_r=-\partial_z A_\phi,\qquad
B_z=\frac1r\partial_r(rA_\phi),\qquad
\nabla\times(\mu_0^{-1}\nabla\times\boldsymbol A)
=J_\phi\boldsymbol e_\phi .
\]

每匝的环向欧姆定律及净电流约束为

\[
\rho J_\phi=\frac{V_n}{2\pi r}-\partial_t A_\phi,
\qquad \int_{A_n}J_\phi\,\mathrm dr\,\mathrm dz=I_n(t).
\]

\(V_n\) 是约束该匝净电流的有效环向电压；接到电源后的实际端电压还包括引线、开关及电路元件。各区由独立电源驱动，区内匝串联，故同一区的匝电流相等（考虑方向符号），而不同区不必相等。

## 2. 磁场、互感和供电电压

磁静态计算先给每匝选择满足 \(\int_{A_n}p_n\,\mathrm dA=1\) 的电流形状 \(J_{\phi,n}=I_np_n\)。固定电阻率下的静态平衡选项为 \(p_n=1/[r h_n\ln(r_{o,n}/r_{i,n})]\)；另一个 `uniform` 选项为 \(1/A_n\)，用于明确规定均匀外加电流密度的比较算例。两者是不同物理源，不应在损耗与场比较时混用。

Biot–Savart 定律给出单个圆环的场核 \(\boldsymbol G_B\)，其中的完整椭圆积分用稳定的互补参数求值。对金属截面求积后叠加所有匝：

\[
\boldsymbol B(r,z,t)=\sum_n I_n(t)
\int_{A_n}p_n(r',z')\boldsymbol G_B(r,z;r',z')\,\mathrm dA'.
\]

同一电流形状定义互感矩阵 \(M_{nm}=\iint p_n p_m M_{\rm loop}\,\mathrm dA_n\mathrm dA_m\) 与储能 \(W_m=\tfrac12\boldsymbol I^\mathsf TM\boldsymbol I\)。代码对自项的可积对数奇点采用 Duffy 变换，对两轴向坐标的积分利用差变量卷积降维，并对近轴和近重合环核作稳定处理。任意电源区映射可归并为区互感矩阵。磁静态结果适于初筛中心场、局部场和低频电压估计。

固定电流形状的集中参数近似为 \(V_k\simeq R_k I_k+\sum_l M_{kl}\dot I_l\)。当皮肤效应与邻近效应使 \(J_\phi\) 随时间重排时，端电压应从时变场与电源约束求得；单个静态 \(M\) 与 \(R\) 不再是完整的端电压模型。因此代码中的电压筛选值不能直接当成已验证的开关或电源额定电压。

## 3. 焦耳热、涡流热和轴对称体积积分

导体局部焦耳热源为 \(q_{\rm ohm}=\boldsymbol J\cdot\boldsymbol E=\rho J_\phi^2\)，单位 \(\mathrm{W/m^3}\)。二维轴对称截面上的**三维总功率**为

\[
P(t)=\sum_n\int_{A_n}2\pi r\,\rho J_\phi^2\,\mathrm dr\mathrm dz,
\qquad Q(t_1,t_2)=\int_{t_1}^{t_2}P(t)\,\mathrm dt .
\]

在 COMSOL 中可选其一：用截面积分并在表达式乘 \(2\pi r\)，或启用“计算体积分/回转几何中的积分”后直接积分 \(q\)。两种约定不能叠加。仓库的独立单环基准将两个积分约定同时显式设置，241 个保存时刻的功率逐点一致。`mf.Qh` 是电磁热源密度，须在相同域、相同回转积分设置下与 \(\rho\,\mathrm{mf.Jphi}^2\) 对照；二者不是两个可相加的热源。[COMSOL 6.4 结果积分设置](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_results.37.027.html)及[电磁损耗变量说明](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_modeling.06.23.html)给出对应软件定义。

### 3.1 可精确分离的传输热

恒定电阻率的单匝直流平衡电流密度、电阻和功率分别为

\[
J_{0,n}=\frac{I_n}{h_n r\ln(r_o/r_i)},\qquad
R_n=\frac{2\pi\rho}{h_n\ln(r_o/r_i)},\qquad
P_{J,n}=R_n I_n^2 .
\]

任意分段线性波形的 \(\int I_n^2\mathrm dt\) 在每段都能解析积分，故固定 \(\rho\) 的传输热不受时间采样积分误差影响。若令 \(J_\phi=J_0+J_e\) 且每匝 \(\int_{A_n}J_e\mathrm dA=0\)，则交叉项

\[
\int_{A_n}2\pi r\,2\rho J_0J_e\,\mathrm dA
=\frac{4\pi\rho I_n}{h_n\ln(r_o/r_i)}
\int_{A_n}J_e\,\mathrm dA=0 .
\]

因此在**同一匝域、同一恒定 \(\rho\)、相同净电流定义**下，\(Q_{\rm total}=Q_J+Q_e\)，其中 \(Q_e=\int\!\int 2\pi r\rho J_e^2\mathrm dA\mathrm dt\ge0\)。这是数学分解，不意味着 COMSOL 自动提供单独“涡流热”变量。若 \(\rho(T,B)\) 改变，需在当前状态重新定义基态及加权正交关系，不能照搬固定 \(\rho\) 的数值拆分。

COMSOL 的 *External Current Density* 特征默认不把外加电流自身的电阻损耗计入电磁损耗；需要另设一致的热源或改用合适的线圈模型。当前 op5 求解模型启用的是八个电流驱动 *Single conductor Coil*，旧的八个外加电流密度节点已禁用；核算时不应重复加 \(I^2R\)。[COMSOL 单导体线圈说明](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.061.html)与[外加电流密度说明](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_electric_fields.07.053.html)可作为软件实现参考。

### 3.2 从缓变近似到有限扩散

把区电流写为向量 \(\boldsymbol I\)，忽略扩散记忆时，投影感应电场给出

\[
P_{e,\rm slow}=\rho^{-1}\dot{\boldsymbol I}^{\mathsf T}
C\dot{\boldsymbol I},\qquad
P_J=\rho\boldsymbol I^\mathsf T G\boldsymbol I .
\]

\(G\) 是各区金属电阻除以 \(\rho\) 的矩阵；\(C\) 是投影后的矢势面积积分，包含不同区之间的交叉耦合，半正定。对于峰值向量 \(\boldsymbol I_*\)、两侧各 \(T\) 的线性斜坡和固定平顶 \(F=1\,\mathrm{ms}\)，有

\[
Q_{\rm slow}(T)=\rho A(F+2T/3)+\frac{2B}{\rho T},\quad
A=\boldsymbol I_*^\mathsf TG\boldsymbol I_*,\quad
B=\boldsymbol I_*^\mathsf TC\boldsymbol I_*,\quad
T_*^{\rm slow}=\sqrt{\frac{3B}{\rho^2A}} .
\]

这说明最优斜坡来自随 \(T\) 增长的传输热与随 \(1/T\) 降低的**缓变**涡流热的平衡。它不是完整磁扩散模型的普适最优时间；当 \(T\) 接近扩散时间，\(1/T\) 渐近式失效，完整能量的最低点可在允许的最短斜坡。以该式批量预筛后，必须用有限扩散模型复算入围波形。

快速的理想薄板模型从 \(\partial_t B=(\rho/\mu_0)\partial_x^2B\) 展开奇数模态 \(q=1,3,\ldots\)，时间常数 \(\tau_q=\mu_0d^2/(\rho\pi^2q^2)\)，求解 \(\tau_q\dot f_q+f_q=\dot B_{\rm app}\)。每个模态在所有斜坡和平顶段解析传播，并将最后时间节点之后的指数衰减热积分到无穷大。矩形两轴薄板之和使用匝中心场，是筛选近似，无法表示完整二维电流重排。原有单极点模型另需几何相关的 FEM 标定。

更高保真度的 `coupled_network` 把每匝截面离散成多个连续环向导电胞元。电磁弱形式为

\[
L\dot{\boldsymbol i}+R\boldsymbol i=C_t^\mathsf T\boldsymbol v,
\qquad C_t\boldsymbol i=\boldsymbol I_t .
\]

\(L\) 是由同一环核积分所得的胞元互感矩阵，\(R\) 是精确环形胞元电阻；\(C_t\) 把胞元电流求和到物理匝。令 \(\boldsymbol i=G_t\boldsymbol I_t+Q\boldsymbol q\)，其中 \(G_t\) 按电阻倒数拆分每匝电流，\(C_tQ=0\)、\(Q^\mathsf TRG_t=0\)。消去未知匝电压得到

\[
(Q^\mathsf TLQ)\dot{\boldsymbol q}+(Q^\mathsf TRQ)\boldsymbol q
=-Q^\mathsf TLG_t\dot{\boldsymbol I}_t .
\]

广义本征模态将其化为独立一阶系统。代码精确传播每个分段线性电流段并积分 \(\rho J^2\)；若脉冲末端所有匝电流已经归零，末端后的剩余感应电流损耗为 \(\tfrac12\sum_m y_m^2(t_{\rm end})\)。近期实现把各模态的时间积分批量化，减少 Python 循环；它不改变物理方程或空间离散精度。对固定几何，先编译互感/本征系统，再扫描电阻率与波形，可避免重复建模。适用范围和四个小环 FEM 对照见[损耗验证记录](loss_theory_validation.md)。

## 4. 从电磁力到结构应力

洛伦兹体力密度为 \(\boldsymbol f=\boldsymbol J\times\boldsymbol B\)。在上述轴对称模型中，\(f_r=J_\phi B_z\)、\(f_z=-J_\phi B_r\)。真正的位移与应力需要对金属、绝缘和约束材料求解动量平衡：

\[
\nabla\cdot\boldsymbol\sigma+\boldsymbol f
=\rho_m\,\partial_t^2\boldsymbol u,
\quad \boldsymbol\sigma=\mathsf C:
[\boldsymbol\varepsilon(\boldsymbol u)-\boldsymbol\alpha_T\Delta T
-\boldsymbol\varepsilon_{\rm pre}] .
\]

准静态近似可在已证明惯性不重要时略去右端。实际计算必须描述层间接触、支撑和轴向端部边界、降温预应力，以及 Zylon/GFRP 的方向相关刚度与强度。\(rJ_\phi B_z\) 是薄环环向载荷量级代理，不是 CuNb 的 von Mises 应力，更不能与 700 MPa 许用值直接比较。当前 Python 包只给力密度/载荷代理；结构强度需独立有限元与材料数据验证。[COMSOL 洛伦兹力说明](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_acdc.26.10.html)提供耦合载荷约定。

## 5. 从热源到温度场与冷却量

瞬态传热方程为

\[
\rho_m c_p(T)\partial_tT
=\nabla\cdot[\mathsf k(T)\nabla T]+q_{\rm ohm}+q_{\rm other} .
\]

金属、环氧、GFRP 和约束材料各有自己的密度、比热和可能各向异性的热导；界面热阻、绝热/对流/液氮边界决定径向及轴向传热路径。除非界面被明确设为绝热且无连接，径向热流不能凭图像颜色认定为零。磁脉冲期间局部热点与整体平均温升不同；\(Q/(mc_p)\) 仅是指定热容、无冷却条件下的**整体平均**估计。液氮蒸发量 \(Q/L_v\) 仅是理想能量账，不能代替瞬时冷却与峰值温度计算。COMSOL 的热源与传热特征可参见[热源说明](https://doc.comsol.com/6.4/doc/com.comsol.help.heat/heat_ug_ht_features.09.016.html)。

完整能量核对至少包括 \(\int P_{\rm source}\mathrm dt\)、各材料内能变化 \(\Delta U=\int_V\int_{T_0}^T\rho_m c_p(\theta)\mathrm d\theta\mathrm dV\)，以及边界流出的热量；若存在其他热源或机械功，也必须列入。Python 包尚无空间传热求解器，温度/应力不能由磁场和总热量自动推出。

## 6. op5 已求解模型对照与未解决项

op5 新几何的历史筛选报告为 340 匝、10 层、8 个独立串联匝电源区、CuNb 金属截面 \(7.5\times12.5\,\mathrm{mm^2}\)、\(18.054+1+18.054\,\mathrm{ms}\) 梯形电流波形。其逐匝解析传输热为 **346.417011 kJ**；旧单极点筛选模型给出总热约 **702.875 kJ**。这两个数字的物理层级不同。

在已求解的 `op5_coupled_v2.mph` 中，对 0–70 ms、281 个时刻的 CuNb 域按回转体积积分，`mf.Qh` 总源热为 **695.621672 kJ**。同一波形的解析传输参考为 **346.416537 kJ**（该时间样本口径），故相减为 **349.205135 kJ**；以求解的 \(J_\phi\) 与基态电流密度逐点作平方拆分，感应偏差平方项 **350.146149 kJ**、交叉项 **−0.941014 kJ**。连续模型的交叉项应为零，这个残差不能直接宣称为真实负涡流热，也不能未经检查就归咎于求解器。需核查每匝电流约束、CuNb 选择域、同一时刻的基态电流及积分与时间采样。

更关键的是，该模型从源热积分得到的 **695.622 kJ** 与 70 ms 的储存热能增量相差 **22.703 kJ（3.26%）**，其来源尚未闭合。旧表中的 **574 kJ** 是没有回转权重的截面积分经时间积分后误标成总焦耳，量纲实际为 J/m；以它减去三维输运热得到的旧 **228 kJ** 也无效。因此约 **1.04%** 的旧解析筛选总热与 FEM 源热差异只能作为**这一几何、这一波形的筛选对照**，不能推广为任意几何、局部温度、应力或材料组合的误差保证。原生模型见 `E:\project\APEX\CuNb-design\APEX_CuNb_B_COMSOL64_168turns_modeling_pack\APEX_CuNb_Optimization\op5_radial_min5\op5_coupled_v2.mph`；仓库的[op5 筛选报告](../studies/apex_cunb_op5/op5_report.md)保留设计参数与筛选假设。

## 7. 复算和验收次序

1. 冻结几何、金属域、每匝净电流、区映射、波形、电阻率和热容的来源与单位；任一项变动都应产生新的算例标识。
2. 确认外加电流密度或线圈特征的实际启用状态。逐匝核对 \(\int_{A_n}J_\phi\mathrm dA=I_n\)，并确认电阻发热只计一次。
3. 明确 `intvolume` 开关，交叉核验 \(\int_V\mathrm{mf.Qh}\mathrm dV=\int_A2\pi r\,\mathrm{mf.Qh}\mathrm dA\)；同域核验 `mf.Qh` 与 \(\rho\,\mathrm{mf.Jphi}^2\)。
4. 分别收敛中心场、局部场、区端压、总源热、涡流增量、脉冲后尾热、局部热点、真实应力；总热接近不能掩盖小分量误差。
5. 对脉冲结束后的全尾部积分，另以源热、储热、边界热流作能量平衡。未闭合时先查定义和耦合，再讨论解析模型误差。

可复算入口：`python -m pytest -q`；单环回转积分审计为 `python validation/comsol/run_validation.py single_a_axisym_audit`。模块细节见[电磁数值推导](first_principles.md)、[损耗有限扩散推导与 FEM 对照](loss_theory_validation.md)和[精度边界](accuracy.md)。

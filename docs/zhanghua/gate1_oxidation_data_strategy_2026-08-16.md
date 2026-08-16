# gate-1（可氧化 Cys）数据策略（2026-08-16）

> 两步模型 stage-1 的训练标签：已知可氧化的 Cys 位点。本文记录当前已建部分、量级现实、
> 扩产路线与可行性判断。配套注册表 `data/registry/gate1_oxidation_sites_v1.tsv`。

## 1. 已建成（Tier-1）：COPLBI 综述 S-亚硝基化列，10 位点 / 6 核验

`data/registry/gate1_oxidation_sites_v1.tsv`，fail-closed 坐标核验：
ACOh4 C172、MEK1 C172、P5CR C5、HA2 C206（番茄）、RAB7 C171、MPK6 C201（拟南芥）
——6 个 verified；GSNOR1 C10、GSNOR C47、LCD C225、PRMT5 C125 为覆盖缺口（蛋白不在
注册蛋白组）。

**诚实量级**：6 个核验位点远不足以训练 stage-1。这是**种子**，不是训练集。stage-1 需要
10²–10⁴ 量级。

### A 步执行结果（2026-08-16，已完成）
- **数据获取通了**：用 Europe PMC `supplementaryFiles` API（绕过 PMC/Wiley 直链拦截）
  下载了番茄 SNO 论文（PMC13205882）与拟南芥 SNO 蛋白组流程论文（PMC12454217）的完整
  补充包（`data/raw/supplements/PXD_tomato_SNO/`）。
- **番茄论文（PMC13205882）补充里没有坐标级 SNO 位点表**——MS 鉴定表只有标准
  Carbamidomethyl/Oxidation，图数据是 GO/DEG 表。此源不产出位点。
- **拟南芥流程论文（PMC12454217）MOESM2 有真 SNO 位点表**：提取 **2974 个拟南芥
  S-nitrosylation 位点**（TAIR 位点 + Cys 位置 + 定位概率），TAIR→UniProt 批量映射
  （2950/2956 位点映射成功），注册蛋白组核验 **2935 通过（98.7%）**。
- **注册表结果**：`gate1_oxidation_sites_v1.tsv` 从 10 行（6 核验）→ **2945 行（2941 核验）**
  ——拟南芥 2939 + 番茄 6，全部 s-nitrosylation。
- 核磺化蛋白组（PMID 39726278）仍非 OA；番茄 SNO 位点表缺位（该论文无坐标级表）。
- 中间产物：`arabi_sno_sites_extracted.tsv` / `tair_to_uniprot_map.tsv` /
  `arabi_sno_sites_verified.tsv`（同目录）。

## 2. 扩产路线（按可行性与机制对齐排序）

### A. 植物 sulfenylation / SNO 文献（机制对齐，量级 10²）
- 植物核磺化蛋白组（PMID 39726278，"The nuclear sulfenome of Arabidopsis"，2025）等
  sulfenylation 数据集——需逐篇提取位点表 + 坐标核验（项目纪律）。
- 拟南芥/水稻 SNO 蛋白组文献（Europe PMC 搜索 412 命中中筛选带位点表的）。
- **判断**：可拿到几百个机制对齐的植物 gate-1 位点；工作量为逐篇提取，属 provenance
  工作流。对 stage-1 仍偏薄，但机制最纯。

### B. UniProt 二硫键注释（量大 10³，但结构性偏置）
- 拟南芥/番茄/水稻参考蛋白组的 DISULFID 特征（UniProt REST，keyword:KW-0182）。
- **判断**：量大但二硫键以分泌蛋白/结构二硫为主，是"永久氧化"而非"调控性氧化"——
  教 stage-1 "胞外 Cys 氧化"不是我们要的 gate-1。可作**明确标注的第二类**，不作主源。
- 技术摩擦：REST `feature:`/`ft_display` 字段语法不稳定（本会话 400 频发），需另走
  JSON 分页提取 features。

### C. 跨物种氧化数据库（量级 10⁴–10⁵）——**降级为显式假设，不作首选路径**
- **dbSNO 已失效**（域名不解析）。替代：
  - 哺乳动物 sulfenylation：Sulfenome / CysOx 等；
  - SNO：PhosphoSitePlus 的 SNO 条目、dbSNO 镜像/论文补充表；
  - UniProt 人/小鼠 DISULFID + SNO；
- **⚠ 与项目注册证据冲突（2026-08-16 修正）**：跨物种迁移在项目已有证据下 **≈0**
  （`multispecies_training_v1_2026-08-10.md` §4.2 LOSO：留一物种全 ≈1.0x base rate；
  解释为"当前特征捕获的是物种内分布记忆，而非通用反应化学"）。**因此"氧化化学跨物种
  保守、人源数据可作 stage-1 先验"是我未经验证的假设，与 LOSO 结论矛盾，不应作为首选。**
- **修正判断**：跨物种 stage-1 是**显式可检验假设**（"gate-1 氧化比 gate-2 持硫化更保守"），
  成立条件是把 stage-1 的判别特征从"物种内分布记忆"升级为通用反应化学（pKa 代理、结构
  反应性）——这正是当前特征没有的。**执行顺序**：先做植物本地 stage-1 基线（A 步），
  证明两步框架本身提效；再用"植物 ± 人源氧化数据"对照检验跨物种 stage-1 是否增值。
  在此之前不引入人/小鼠蛋白组基础设施（=新预注册 + 大量投入，而假设未被支持）。

## 3. 推荐顺序（务实，2026-08-16 修正后）

1. **A 步（首选，唯一有项目证据支持的方向）**：挖 1–2 篇**植物** sulfenylation/SNO
   数据集（核磺化蛋白组 PMID 39726278 + 一篇拟南芥 SNO），把 Tier-1 扩到 ~100–300 个
   核验位点。依据："同物种数据是唯一实证有效方向"（LOSO 结论）——gate-1 先做植物本地。
2. **B 步（辅助类）**：UniProt 二硫键作明确标注的第二类，不作主源（结构二硫 ≠ 调控氧化）。
3. **C 步（显式假设，非首选）**：跨物种 stage-1 需先证明两点——(i) 两步框架在植物本地
   基线上确实提效；(ii) 人源氧化数据在"植物 ± 人源"对照中增值。满足前才立项（新预注册 +
   人/小鼠蛋白组基础设施）。
4. 每一步都在 stage-1 模型上验证是否提升蛋白内 Hit@2（金标准对照），不盲目堆数据。

## 4. 与两步模型的连接

`P(持硫化) = P(可氧化) × P(持硫化) / P(可氧化)`（因持硫化 ⊆ 可氧化）——stage-1 的
价值是**数据增强**（用 10²–10⁵ 氧化数据把"哪些 Cys 反应性"学好），不是分解本身。
当前 Tier-1（6 位点）无法支撑该价值；A 步 + C 步之后才可检验"两步模型真的提效"。

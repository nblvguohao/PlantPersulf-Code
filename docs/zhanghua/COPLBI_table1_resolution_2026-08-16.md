# COPLBI-D-26-00068 Table 1 位点逐条解析（2026-08-16）

> 目的：把综述 Table 1 两列（持硫化 / S-亚硝基化）的每个位点按**原始论文来源**解析到注册
> 参考蛋白组的正确 accession，并核验该位置的 Cys。第一轮审计曾错误地对全部位点假设拟南芥
> 蛋白组——综述 Table 1 实际跨物种（番茄 / 拟南芥 / 水稻），这是"9/10 SNO 核验失败"的
> 根因。本文记录每条解析结果与来源。审计输出：
> `scripts/evaluate_oxiptm_sites.py` → `results/diagnostics/oxiptm_sites_v1.json`。

## 持硫化列（10 基因）

| 基因/位点 | 物种（综述正文） | 解析 accession | Cys@位点 | 状态 |
|---|---|---|---|---|
| DES1 C44/C205 | 拟南芥 | F4K5T2 | ✓/✓ | **verified**（Shen 2020, Plant Cell） |
| RBOHD C825/C890 | 拟南芥 | Q9FIJ0 | ✓/✓ | **verified**（Shen 2020） |
| SnRK2.6 C131/C137 | 拟南芥 | Q940H6 | ✓/✓ | **verified**（Chen 2020/2021, Mol Plant） |
| ATG4a C170 | 拟南芥 | Q8S929 | ✓ | **verified**（Laureano-Marin 2020, Plant Cell） |
| ABI4 C250 | 拟南芥 | A0MES8 | ✓ | **verified**（Zhou 2021, Mol Plant） |
| PAD3 C440 | 拟南芥 | Q9LW27 | ✓ | **verified**（Zhang 2026, PCE） |
| CAT1 C234 | 番茄 | P30264 | ✓ | **verified**（Li 2020, PPB） |
| bZIP68 C171 | **水稻**（Ma X 2026, IJMS 27:3841, 开放获取全文核验） | **A2YXP7**（OsI_30119, 435aa, chr8） | C@171 + **C@245** 双 Cys 与论文两处修饰位点完全一致 | **verified**（=NCBI 125562410；注册蛋白组中唯一具此双 Cys 模式的 BZIP；注：条目为 indica 注释，japonica 位点 ID 未解析出——Ensembl 503、NCBI 无 gene 记录——但为参考中该 bZIP68 位点的唯一代表） |
| APX1 C168 | 番茄 | — | — | **unmappable**（known_controls 已记录，番茄候选位 168=Ala） |
| POD5 C61 | 番茄 | — | — | **unmappable**（番茄零命中） |

## S-亚硝基化列（10 基因）

| 基因/位点 | 物种（综述正文） | 解析 accession | Cys@位点 | 状态 |
|---|---|---|---|---|
| ACOh4 C172 | **番茄**（Liu 2023, New Phytol 239:159-173） | **A0A3Q7FZA2**（366aa） | C@172 ✓ | **verified**：NCBI Gene **LOC101265426** = "1-aminocyclopropane-1-carboxylate oxidase homolog 4"（蛋白 460382410）与参考序列**精确一致**（参考头注为"Fe2OG dioxygenase"，名称扫描漏掉的原因） |
| MEK1 C172 | **番茄**（Fang 2026） | O48616（+A0A3Q7JS13） | C@172 ✓ | **verified**（SlMEK1, MAPKK） |
| GSNOR1 C10 | 拟南芥（Zhan 2018） | —（仅 Q0WM36 195aa 片段 P@10） | ✗ | **coverage gap**：蛋白组 v2 缺规范 GSNOR1 |
| GSNOR C47 | **番茄**（Huang 2026） | — | — | **coverage gap**：番茄 GSNOR 不在蛋白组 |
| LCD C225 | **番茄**（Huang 2026） | — | — | **coverage gap**：番茄 L-CD 不在（仅 D-CD） |
| P5CR C5 | **番茄**（Liu 2024） | A0A3Q7FME1（SlP5CR） | C@5 ✓ | **verified**（SlP5CRC5S 突变系） |
| HA2 C206 | **番茄**（Wei 2025, Plant Cell 37:koaf035） | **Q9SPD5**（GN=LHA2, 956aa） | C@206 ✓ | **verified**：与 NCBI "plasma membrane H+-ATPase isoform LHA2"（5901757）**100% 一致**；704aa 旧截短记录 P23980（G@206）为同基因的过期条目 |
| RAB7 C171 | 拟南芥（Lin 2023） | Q9XI98（RABG3E） | C@171 ✓ | **verified**（首轮误用 RAB7A A0A1P8AXJ0） |
| PRMT5 C125 | 拟南芥（Hu 2017） | — | — | **coverage gap**：PRMT5/SKB1 不在蛋白组 v2 |
| MPK6 C201 | 拟南芥（Wang 2025） | Q39026 | C@201 ✓ | **verified** |

**关键更正**：首轮 9/10 SNO "失败"中，6 个实际是番茄位点（MEK1、P5CR、GSNOR、LCD、HA2、
ACOh4），我误在拟南芥蛋白组里查。逐条 primary-source 解析后 **SNO 6/10、持硫化 11/13
可验证**；剩余不可验证者全部是注册蛋白组的**覆盖缺口**（GSNOR1 仅 195aa 片段、番茄 GSNOR、
番茄 LCD、PRMT5/SKB1 不在蛋白组），不再是位点核验失败。APX1/POD5 番茄不可解（已知）。

## 判别测试（解析后，6 SNO + 8 持硫化蛋白）

- 蛋白内排名（per-site）：持硫化 `[3,3,2,4,6,3,2,7,1,3,8]` 均值 **3.82**（11 位点 / 8 蛋白）；
  SNO `[4,4,1,1,1,1]` 均值 **2.0**（6 位点 / 6 蛋白；P5CR 为 1-Cys 蛋白 rank 1 平凡值，
  排除后 SNO 均值 2.2）。
- 观测均值差（持硫化 − SNO）= **+1.82**；单侧置换（持硫化更低/更优）**p=0.976**（B=999, seed 20260816）。
- **读数：冻结持硫化训练模型在蛋白内系统性不把持硫化位点排到 S-亚硝基化位点之上——方向
  明确相反**（ACOh4/HA2/MEK1 三例 SNO 位点均蛋白内 rank 1，持硫化组无一 rank 1 之外优势）。
  在两步框架下与"模型主要捕获 gate-1（通用可氧化 Cys）而非 gate-2（持硫化特异）信号"一致；
  也是共肽零结果的平行证据。小样本，描述性，非注册端点。

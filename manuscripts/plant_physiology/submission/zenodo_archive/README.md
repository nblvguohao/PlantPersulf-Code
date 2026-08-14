# Zenodo 数据与代码归档 — 操作指南

目标：把论文承诺「公开归档」的全部工件发布到 Zenodo，获得永久 DOI。
归档内容与手稿 Data availability 一节一一对应。

## 1. 归档内容清单

| 内容 | 路径 | 大小 | 说明 |
|---|---|---|---|
| 冻结发布包（完整） | `results/candidates/multispecies_v2_candidate_release_v1/` | ~8.5 MB | 模型权重、候选表 179,736 行、匹配对照、签署协议、功效表、SHA256SUMS |
| 数据归档 | `results/data_archive/nc_entry_gates/` | ~68 KB | 排除台账（含 EX-009）、输入清单、负面结果登记（含 pLDDT 纠正）、复现命令与环境锁 |
| 保守性/结构分析输出 | `results/cross_species_conservation/` | ~56 KB | `conservation_v2.json`、`structural_context_v2.json`（含 v1 历史版） |
| 补充表 S1 | `manuscripts/plant_physiology/supplements/` | — | `supplemental_table_s1.pdf` + `.tsv` + `.tex` |
| 代码快照 | 见下方命令 | — | 冻结标签 `gate0-release-v1-20260813` 处的仓库快照 |

打包命令（PowerShell，仓库根目录）：

```powershell
$dst = "manuscripts/plant_physiology/submission/zenodo_archive/files"
New-Item -ItemType Directory -Force $dst | Out-Null
Copy-Item -Recurse "results/candidates/multispecies_v2_candidate_release_v1" "$dst/release_bundle" 
Copy-Item -Recurse "results/data_archive/nc_entry_gates" "$dst/data_archive"
Copy-Item -Recurse "results/cross_species_conservation" "$dst/cross_species_conservation"
Copy-Item "manuscripts/plant_physiology/supplements/supplemental_table_s1.pdf","manuscripts/plant_physiology/supplements/supplemental_table_s1.tsv" $dst
git archive -o "$dst/code_snapshot_gate0.zip" gate0-release-v1-20260813
Compress-Archive -Path "$dst/*" -DestinationPath "$dst/../zenodo_archive_v1.zip" -Force
```

（`files/` 与 zip 为提交物副本，勿入库；README 本文件入库。）

## 2. Zenodo 操作步骤

1. 登录 <https://zenodo.org>（PI 账号；如无账号先注册）。
2. **New upload** → 拖入 `zenodo_archive_v1.zip`（Zenodo 会在 Publish 时解包为文件树）。
3. 元数据填写：

| 字段 | 建议值 |
|---|---|
| Title | `PlantPersulf: multispecies cysteine persulfidation dataset, frozen tomato candidate release, and cross-kingdom conservation analyses (release v1)` |
| Creators | 手稿作者表 7 人，顺序一致 |
| Description | 手稿 Data availability 一节内容 + 本 README 第 1 节清单 |
| License | **需要你拍板**：数据/文档建议 CC-BY-4.0；代码建议 MIT（仓库当前无 LICENSE 文件，建议在发布前补上，见 4） |
| Keywords | persulfidation; hydrogen sulfide; cysteine modification; tomato; preregistration; PU learning |
| Related identifiers | OSF 预注册 DOI（选 `isSupplementTo`/`isDerivedFrom` 关系，指向最终论文 DOI 可后补） |
| Version | v1 |

4. **Publish** → 获得 DOI（形如 `10.5281/zenodo.XXXXXXX`）。注意：Publish 后不能改文件，修改需走 New Version（会生成新 DOI）——所以发布前先核对 SHA256。

## 3. 发布后动作

- 把 Zenodo DOI 发给我，填入稿件的 Data availability 两处「Zenodo」位置并重建 PDF。
- 同步填写 `../README.md` 总清单。

## 4. 发布前检查

- [ ] 仓库补 LICENSE（MIT for code；CC-BY-4.0 声明放在 Zenodo 元数据即可）
- [ ] 排除台账中是否还有未了结条目（`results/data_archive/nc_entry_gates/exclusions/excluded_samples.tsv` 应与 EX-001…EX-009 对应齐全）
- [ ] `environment_lock.txt` 与 `commands.md` 已随包上传
- [ ] 代码快照确为 `gate0-release-v1-20260813` 标签（`git archive` 自动保证）

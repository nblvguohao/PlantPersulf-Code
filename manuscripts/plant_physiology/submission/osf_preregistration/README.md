# OSF 预注册 — 操作指南

目标：把已共同签署的 SAP 协议与统计分析计划登记到 OSF，获得预注册 DOI。
登记内容与冻结包完全一致（`sap_protocol.json` / `analysis_plan.json`，状态
`signed_2026-08-13`），登记动作**不改动**任何冻结产物。

## 1. 上传文件清单（全部来自冻结发布包）

| 文件 | 路径 | 说明 |
|---|---|---|
| `sap_protocol.json` | `results/candidates/multispecies_v2_candidate_release_v1/` | SAP 协议 v1.0（签署版） |
| `analysis_plan.json` | 同上 | 统计分析计划 v1.0（签署版） |
| `co_signature_package.md` | 同上 | 人类可读签署摘要 |
| `power_table.json` | 同上 | 144 格功效网格（分析计划引用） |
| `SHA256SUMS` | 同上 | 全部冻结产物哈希 |

打包命令（PowerShell，仓库根目录）：

```powershell
$src = "results/candidates/multispecies_v2_candidate_release_v1"
$dst = "manuscripts/plant_physiology/submission/osf_preregistration/files"
New-Item -ItemType Directory -Force $dst | Out-Null
Copy-Item "$src/sap_protocol.json","$src/analysis_plan.json","$src/co_signature_package.md","$src/power_table.json","$src/SHA256SUMS" $dst
Compress-Archive -Path "$dst/*" -DestinationPath "$dst/../osf_preregistration_v1.zip" -Force
```

（`files/` 与 zip 均为提交物副本，勿入库；README 本文件入库。）

## 2. OSF 操作步骤

1. 登录 <https://osf.io>（PI 账号，建议项目负责人持有）。
2. **新建项目**，建议标题：
   `Site-Aware Prioritization (SAP) of Tomato Cysteine Persulfidation Sites: Preregistered Blind-Cohort Validation (v1, signed 2026-08-13)`
3. 在项目 **Files** 页上传第 1 节全部 5 个文件。
4. 点击 **Registrations → New Registration**，表单选 **OSF Preregistration** 模板。
5. 表单各栏填写（内容全部摘自上传文件，不要改写数字）：

| OSF 表单栏 | 填写内容 |
|---|---|
| Study Information → Title | 同上标题 |
| Study Information → Description | `co_signature_package.md` 第 1 节「冻结了什么」摘要 |
| Hypotheses | SAP `primary_endpoint` 原文：Top-200 候选相对匹配对照的确认富集（单侧 Fisher 精确检验，α=0.05，OR>1，≥2 个确认候选位点） |
| Study design | SAP `blind_design`：593 位点盲法检测、盲标识列表、哈希锁定的标识-标签映射 |
| Analysis plan | `analysis_plan.json` 的 `primary_analysis` / `secondary_analyses` / `unblinding_procedure` 原文 |
| Other | 注明：「Complete signed protocol and analysis plan attached as JSON (v1.0); these files take precedence over this form.」 |

6. **公开性选项**：预注册登记通常立即公开。如担心被抢先，可选 **Embargoed registration**（OSF 支持最长 4 年禁运），并用「View-only link」提供给编辑部/审稿人核验——PP 只需要 DOI 与可核验性，两种方式都满足。
7. 点击 **Register**（不可逆）。记录 DOI，形如 `10.17605/OSF.IO/XXXXX`。

## 3. 注册后动作

- 把 DOI 发给我，我将其填入手稿的 5 处「preregistration DOI to be inserted」位置（摘要注释、R5、R7、Methods、Data availability）并重建 PDF。
- 同步填写 `../README.md` 总清单中的对应勾选项。

## 4. 边界提醒

- 登记版文件必须与冻结发布包字节一致（`SHA256SUMS` 已含）。
- 若日后需要修订协议：按签署机制走「有日期、带哈希的重新签署」，OSF 上登记新版本，**不修改**已登记版本。

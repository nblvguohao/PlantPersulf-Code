# 竞争工具数据可获取性审计（2026-08-10）

> 目的:评估"与 pCysMod / Sul-BertGRU 互跑公平比较"的可行性(用户决策:
> 本阶段只审计数据,不跑模型)。审计对象:数据可下载性、数据与论文的
> 一致性、物种构成、与 Seville 基准的重叠。
> 方法:站点直接抓取(证书过期,curl -k 绕过)、GitHub API、UniProt
> REST(accession→物种/序列映射)。所有原始文件按需重取,中间产物不入库。

## 1. Sul-BertGRU(Bioinformatics 2025, btaf078)——数据可获取 ✅,但有三处不一致

数据源:`https://github.com/Severus0902/Sul-BertGRU`(论文 Data availability
声明),文件 `data/positive1.txt`(108 KB)、`data/negative1.txt`(668 KB),
FASTA 格式(>accession + 31 残基肽段,窗口外补 `X`)。

### 1.1 行数与论文一致,但"中心 Cys"只有 1,702/2,705

| 检查项 | 论文声明 | 下载文件实测 | 判定 |
|---|---|---|---|
| 阳性记录数 | 2,705 | 2,705(5,410 行 FASTA) | ✅ 一致 |
| 阴性记录数 | 16,697 | 16,697(33,394 行) | ✅ 一致 |
| 窗口 | 31 aa,Cys ±15 | 长度 12–31 不一;仅 1,702/2,705 记录中心位为 C | ⚠️ **不一致** |
| 负样本构造 | 论文:16,697 全量 | README 却写"随机选与阳性等量的阴性"(与论文矛盾) | ⚠️ README 与论文矛盾 |

结论:**其 2,705 阳性集的构造在下载文件中无法以"中心 Cys 位点"标准
重现**(1,003 条记录中心非 C;序列未严格居中/对称)。若互跑,必须先用其
`spilt_seq.py` 复现预处理,否则无法确知训练时用的窗口。

### 1.2 物种构成:78% 人类 + 21% 拟南芥(全量 1,702 唯一 accession 映射)

| 物种 | 蛋白数 | 占比 |
|---|---|---|
| Homo sapiens | 1,335 | 78.4% |
| Arabidopsis thaliana | **364** | 21.4% |
| 未映射 | 3 | 0.2% |

正文未报告物种;实测**超过 1/5 是拟南芥**——植物位点并不像正文暗示的
那样缺席,且拟南芥 persulfidation 的公开位点来源正是 Seville 体系,
同源风险高(见 1.3)。负样本抽样(300):84% 人类 + 16% 拟南芥。

### 1.3 与 Seville 基准的重叠:蛋白级 33,位点级 0(有审计局限)

- 将 364 个拟南芥蛋白的阳性肽段映射回 UniProt 全长序列(去 X 子串匹配):
  **569 个位点映射成功**;
- 与我们基准(350 阳性蛋白)比对:**蛋白级交集 33 个,位点级交集 0 个**;
- **审计局限**:1,338 个肽段未能匹配当前 UniProt 序列(大概率 iCysMod
  收录时使用了旧版序列)——无法排除未匹配肽段中存在与基准重叠的位点,
  位点级 0 重叠是"至少 569 个映射位点无重叠"的保守下限,不是全量结论。
  33 个蛋白共享但位点不同,亦不能解释为同源位点。

## 2. iCysMod(数据源本体)——表单门控,未验证实际交付

- 站点:icysmod.omicsbio.info,证书过期但可访问;数据库论文
  Wang et al., *Brief Bioinform* 2021, doi:10.1093/bib/bbaa400
  (85,747 位点 / 31,483 蛋白 / 48 真核物种 / 8 类 Cys 修饰,
  S-sulfhydration 为其中一类);
- Download 页需**提交个人信息表单**获取链接,页面注明数据集在论文发表后
  可下载;未实际提交表单验证交付(隐私与成本考量,记录为未验证项);
- 附带发现:iCysMod 摘要报告 37,841 个事件涉及 119 种"同一 Cys 位点
  PTM 共现"——与我们的 crosstalk grammar 主题直接相关,值得后续引用。

## 3. pCysMod(Front Cell Dev Biol 2021, 617366)——训练数据当前不可获取 ❌

- 论文声称"proteins and peptides used in this study were uploaded in the
  web server and users can download the relevant data in the 'Help' section";
- 实测:站点导航中 Download 链接已**注释掉**;download.php 页面为占位壳
  ("Loading... 0%",无任何链接);Help 页仅含预测使用说明;
- 结论:pCysMod 的 5 类修饰训练数据**当前无法下载**——互跑 pCysMod
  在数据侧不可行(最多在其 web server 上做定性预测对比,拿不到训练集)。

## 4. 互跑可行性结论

| 维度 | Sul-BertGRU | pCysMod |
|---|---|---|
| 数据获取 | ✅ GitHub 直下 | ❌ 下载入口失效 |
| 数据-论文一致性 | ⚠️ 中心 C 1,702/2,705;README 矛盾 | — |
| 环境可复现性 | ⚠️ Python 3.6 / TF 1.15 / PyTorch 1.10,无 requirements.txt | 未提供代码(仅 web server) |
| 与基准重叠风险 | ⚠️ 拟南芥子集 21%,蛋白级重叠 33 | 数据不可得,无法评估 |

**建议**:互跑仅对 Sul-BertGRU 有意义且需满足三个前置:
(1) 用其 `spilt_seq.py` 复现 2,705 阳性集的构造并核对中心 C 差异;
(2) 按物种分层报告(人类/拟南芥),拟南芥子集与 Seville 数据同源风险
    使跨物种结论不可靠;
(3) 环境隔离(老版本 TF/PyTorch 与当前项目冲突)。鉴于当前阶段目标是
    literature-comparable 数字(新轨道 `pu_ranker_protein_split_v1` 已
    落地,见 two-track policy 更新),互跑推迟至有张华团队独立队列后
    一并评估。

## 5. 对 phase_z 审计的影响

本审计实证了 `docs/phase_z_evidence_audit.md` 的既有结论:"非 Seville
植物位点级数据"全领域稀缺——连 Sul-BertGRU 的拟南芥子集都无可证实的
位点级重叠,且 iCysMod 收录的拟南芥序列版本与当前 UniProt 差异明显
(1,338/1,907 肽段无法映射),进一步说明该领域位点数据版本混乱、难以
交叉核验。这强化了与张华团队合作中"位点级数据 + 完整 provenance"的
稀缺价值。

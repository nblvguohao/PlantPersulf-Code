# Phase A1 — Persulfidation Site-Evidence Source Map & Acquisition Design

> **Cycle role:** This is Phase A1 of the post-Gate-1 roadmap
> (`docs/superpowers/plans/2026-07-21-post-gate1-roadmap.md`). It is a **source-map
> and acquisition design only**. It writes no parser, creates no benchmark label,
> and upgrades no evidence class. Its single job is to decide *where* real
> Cys-site persulfidation evidence for each discovery/external study actually
> lives, at what data level, whether it is legally obtainable, and how to register
> it — so that Phase A2 can parse it under TDD and Phase A3 can re-run the
> readiness gate.

---

## 中文执行摘要

- **背景**：Task 5 readiness gate 判定 **STOP**（0 个合格位点、0 个合格研究，需 ≥2）。根因是四个硫巯基化研究目前只解析了外围文件，**位点级 persulfidated-Cys 证据尚未触达**。
- **本文档任务**：为 PXD006140 / PXD024061 / PXD035795 / PXD039999 各画一张"位点证据在哪"的地图，标注数据等级、获取路径、开放许可、以及**能否真正给出位点级（而非仅蛋白级）证据**这一决定性未决问题。
- **三条获取路线（按可行性排序）**：
  1. **Route S（推荐首选）**：同行评议**补充表**（Data 等级 C）——四篇论文全部有 DOI 且期刊均开放获取，补充表是最干净的位点来源，无需重跑 GB 级原始数据。
  2. **Route M（唯一已在手）**：PXD035795 已下载的 `mzid`（C 残基已定位修饰）+ 已登记的方法 PDF → 只差一条经审阅的"质量/化学标签 → persulfidation 位点"映射规则。
  3. **Route H（重、暂缓）**：解析 GB 级搜索输出（PXD024061 的 2.42GB MaxQuant `txt` zip、PXD006140 的 Mascot `.dat`）——权威但工程量大，仅作交叉验证。
- **决定性风险**：若补充表只给蛋白级富集、不给具体 Cys 位点，则 Route S 也无法产出 `site_ms`，项目按契约走 **Phase Z**（降级为数据资源/证据审计）。**这一点必须由 A2 用真实文件实证，不得假设。**
- **达成 ≥2 研究的最可能组合**：PXD035795（Route M）+ {PXD006140 | PXD024061 | PXD039999} 之一的补充表（Route S）。

> **⚠️ 2026-07-21 实测结论（见 §5「实测日志」，含 §5.5 / §5.6 深挖）**：已实际下载并检视四篇论文的公开产出 + PXD035795 mzid/方法 PDF + PXD024061 的 MaxQuant 搜索输出。**两个研究均确认有真实位点级 persulfidation 证据 → Gate 1 可从公开数据达成 GO,Phase Z 不再是默认路线。**
> **① 研究 #1 = PXD006140（§5.5）**：Dataset S3 用作者自定义的 Cys persulfidation 修饰(`MOD:99998 Sulfide` / `99997 CN-Biotin-Sulfide`),FDR<1%,约 **356 个位点**(323 单蛋白无歧义 + 51 isoform 多映射入冲突表);坐标经真实 UniProt 交叉验证(O03042 C284、Q944G9 C79)。
> **② 研究 #2 = PXD024061（§5.6）**：其 `txt_persulfproject.zip` 里的 MaxQuant `Sulfide(C)Sites.txt` / `CianoBiotin(C)Sites.txt` 是标准位点表,去 decoy/contaminant、按定位概率 ≥0.75 得 **~76 个 class-I 位点**(全部 Amino acid=C)。仅用 HTTP Range 取了这两个小成员,未下 2.42GB 整包。
> **③ 其余**：PXD039999 公开材料纯蛋白级;PXD035795 的 MS 标签(DCP=磺烯酸反应性、NBF=通用巯基)非作者认定的 persulfidation 位点标记,且其 New Phytol 文章**不在 PMC/非 OA**(需机构权限)——但 #1+#2 已够 Gate 1,#2 不再依赖它。
> **④ 重要限制**：两个研究同属 Seville(Romero/Gotor)实验室、同套 tag-switch 化学 → 化学方法**不独立**,leave-study-out 可做但方法多样性有限,须在结论中显式声明(见 §5.6)。
> **A2 首个目标 = PXD006140 Dataset S3 位点解析;第二目标 = PXD024061 MaxQuant Sites 位点解析。**

---

## 1. Scope and non-goals

**In scope (A1):**
- Enumerate every candidate site-evidence carrier per study, grounded in the
  registered `data/registry/files.tsv` inventory and the official publication DOI.
- Classify each candidate by expected data level (A–E per the Codex contract) and
  acquisition tractability.
- Specify the provenance-registration schema for any newly acquired supplementary
  source, reusing the `evidence_methods.tsv` precedent already set for PXD035795.
- State, per study, the **decisive open question**: does the cleanest obtainable
  source carry residue-resolved Cys evidence, or only protein-level enrichment?

**Explicit non-goals (deferred to A2/A3):**
- No parser, no `site_ms`/`site_mutagenesis`/`site_biochemical` assignment.
- No download of multi-GB raw/search files into the scientific pipeline.
- No benchmark, label, split, feature, or model artifact.
- No claim that any candidate *is* a positive until A2 verifies it against the
  UniProt sequence and (for MS evidence) a registered method mapping.

---

## 2. Study publications (verified from registered PRIDE metadata)

All four studies are from the Romero/Gotor group (CSIC–Universidad de Sevilla).
DOIs and PubMed IDs below were read from the registered
`data/registry/cache/pride/<ACC>.json` `project.references` field (not guessed).

| Accession | Publication | DOI | PMID | Journal / OA status |
|---|---|---|---|---|
| PXD006140 | Aroca, Benito, Gotor, Romero 2017 | `10.1093/jxb/erx294` | 28992305 | J Exp Bot (OA supplement) |
| PXD024061 | Jurado-Flores, Romero, Gotor 2021 | `10.3390/antiox10040508` | 33805243 | Antioxidants (MDPI, CC-BY) |
| PXD035795 | García-Calderón et al. 2023 | `10.1111/nph.18838` | — | New Phytologist (supplement) |
| PXD039999 | Jurado-Flores, Aroca, Romero, Gotor 2023 | `10.1093/jxb/erad165` | 37148339 | J Exp Bot (OA supplement) |

**License note:** PXD024061 (MDPI CC-BY) is unambiguously redistributable.
J Exp Bot and New Phytologist supplementary files are publicly downloadable; the
registry must record each file's specific license/usage string as observed at
download time (do not assume CC-BY for non-MDPI supplements).

---

## 3. Per-study site-evidence source map

Legend for **Site level?** column: **YES** = residue-resolved Cys localization
expected; **MAYBE** = must be verified by inspection in A2; **NO** = protein- or
identification-level only.

### 3.1 PXD006140 — discovery persulfidation (Arabidopsis, tag-switch/TMT)

Registered repository files (`files.tsv`):

| File | Category | Size | Site level? | Notes |
|---|---|---|---|---|
| `omssa.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt` | OTHER | 18 MB | NO | Already downloaded; parser yields only `psm_coordinate_only`; the 3 Cys rows are Met-oxidation-adjacent coordinates, **not** persulfidation calls |
| `omssa.ne.…mgf.txt` | OTHER | 18 MB | NO | Same as above |
| `F105505.AAroca.mods.dat` | SEARCH | 221 MB | MAYBE (heavy) | Mascot `.dat` with modification search; residue localization present but requires a Mascot-`.dat` parser (Route H) |
| `F105506.AAroca.dat` | SEARCH | 219 MB | MAYBE (heavy) | Mascot `.dat` |
| `…mgf.K/T.…t.xml` (×4) | OTHER | 220–370 MB | MAYBE (heavy) | TPP/pepXML-like; heavy |
| `…cmpd.mzXML` | RAW | 199 MB | NO | Peak list |

- **Primary route (S):** Aroca 2017 J Exp Bot supplementary tables. The paper's
  comparative persulfidation proteome reports persulfidated proteins/peptides;
  **A2 must confirm whether the supplement lists specific persulfidated Cys
  positions or only proteins.** This is the decisive open question for this study.
- **Secondary route (H, deferred):** parse `F105505.AAroca.mods.dat` (Mascot) for
  residue-localized persulfidation-tag modifications. Authoritative but GB-scale
  and parser-intensive; use only as cross-check if the supplement is site-level.
- **Reference sequences already registered:** Q93VK9, Q9ZW96 (coordinate-validation).

### 3.2 PXD024061 — external, nitrogen starvation (Arabidopsis, MaxQuant, label-free)

| File | Category | Size | Site level? | Notes |
|---|---|---|---|---|
| `checksum.txt` | OTHER | 3 KB | NO | Downloaded; provenance inventory only |
| `txt_persulfproject.zip` | SEARCH | **2.42 GB** | MAYBE (heavy) | MaxQuant `txt` output; may contain `modificationSpecificPeptides.txt` / a Cys-site table (Route H) |
| `arabidopsis_uniprot_072020_identified.fasta` | FASTA | 4.6 MB | — | Search DB; useful for ID/version mapping |
| `*.raw` (×24) | RAW | ~2.5 GB each | NO | Thermo raw |

- **Primary route (S):** Jurado-Flores 2021 *Antioxidants* supplementary tables
  (**CC-BY, cleanest license**). MDPI supplements are directly downloadable.
  A2 must confirm site-vs-protein granularity.
- **Secondary route (H, deferred):** stream only the site-relevant members of
  `txt_persulfproject.zip` (e.g. a MaxQuant sites/modificationSpecificPeptides
  table) without materializing the full 2.42 GB into the pipeline; register the
  extracted member's own SHA256.

### 3.3 PXD035795 — external, non-photorespiratory (Arabidopsis, dimedone/DCP + NBF) — **CLOSEST TO SITE-LEVEL**

| File | Category | Size | Site level? | Notes |
|---|---|---|---|---|
| `peptides_1_1_0.mzid.gz` | RESULT | 1.9 MB | **MAYBE (in hand)** | Downloaded & audited: 124 Modification elements, residues blank/**C**/K, masses 163.0012 / 168.0786 / 196.08 / 394.1557, CV term only `MS:1001460 unknown modification` |
| `peptide.csv` | OTHER | 1.7 MB | NO | Identification-level; free-text PTM summary, no verified residue location |
| `proteins.csv` | OTHER | 0.5 MB | NO | Protein-level only |
| `SDRF.txt` | EXP DESIGN | 5 KB | — | Names search mods NBF_N/NBF_K/NBF_C/DCP/dcp-ac |
| `*.mgf` / `*.raw` / `*.mzxml` | PEAK/RAW | 57 MB–772 MB | NO | Spectra |

- **Primary route (M, already in hand):** the registered mzIdentML localizes
  modifications on **C** residues. The one missing link is a **reviewed method
  mapping** that ties a specific mass/chemical tag (e.g. a DCP-derived delta on C)
  to persulfidation-site evidence. The method source (`Persulfidation_protects.pdf`,
  DOI 10.1111/nph.18838) is **already registered** in `evidence_methods.tsv`.
- **Evidence-policy gate (unchanged from content-audit design):** `site_ms`
  requires (1) a threshold-passing mzid identification, (2) exact peptide/protein/
  location/Cys/SHA256, (3) a registered method source explicitly mapping the
  observed tag/mass to persulfidation, (4) no unresolved conflict, (5) a
  reproducible locator. Absent (3), records stay `unresolved`/`identification_only`.
- **Supplementary route (S) also available:** New Phytol supplement for a curated
  site table if the mzid mapping proves ambiguous.

### 3.4 PXD039999 — external, drought + H₂S (Arabidopsis, timsTOF `.d`)

| File | Category | Size | Site level? | Notes |
|---|---|---|---|---|
| `checksum.txt` | RAW | 1.8 KB | NO | Downloaded; provenance only |
| `proteinSeq.txt` | SEARCH | 14 MB | — | Search DB sequences |
| `*.d.7z` (×6) | RAW | ~8 GB each | NO | Bruker timsTOF raw |
| `*_5.3.556.mgf` (×6) | PEAK | ~0.5 GB each | NO | Peak lists |

- No processed site/result table exists in the repository (Stage A confirmed empty
  of results). **Only viable near-term route is S:** Jurado-Flores 2023 J Exp Bot
  supplementary tables. A2 must confirm granularity.
- Route H (reprocessing 48 GB of Bruker raw) is **out of scope** for the MVP.

---

## 3.5 Blocker taxonomy — *why* Gate 1 is STOP (it is not a parsing bug)

A recurring misframing is that the STOP is a "dataset parsing mismatch" to be
debugged. It is not. The parsers behave correctly on what they are given; the
blocker is a layered **evidence-availability** problem. Four distinct blocker
types, none of which is a parser defect:

| Type | Name | Meaning | Fixable by… | Example |
|---|---|---|---|---|
| **A** | Coverage gap | The site-bearing file was simply never downloaded/parsed | Acquisition | PXD024061/PXD039999 had only `checksum.txt` in Stage A |
| **B** | Source selection | A file *was* parsed, but it is the wrong carrier for site evidence | Choosing the right source | PXD006140 OMSSA `.mgf.txt` (spectra) vs the search `.dat` / supplement |
| **C** | Scientific classification (fail-closed) | Data present, residue localization present, but the modification is **not author-designated as a persulfidation site** | Reviewed method judgment | PXD035795 mzid tags DCP/NBF (see §3.3) |
| **D** | Public-data ceiling | The residue-resolved persulfidation **site** is *not published* in any public output of the study | **Nothing** — inherent limit; triggers Phase Z | Curated persulfidation outputs are protein-level (see §5) |

Every gate the pipeline enforces — refusing to call Met-oxidation a
persulfidation, keeping unlocalizable mods `unresolved`, failing closed on
missing files — is the integrity contract **working as designed**, not a bug.
Type **D** is the decisive one: it cannot be engineered away and is the direct
trigger for the Phase Z pivot.

---

## 4. Acquisition-route ranking (project-wide)

| Route | What | Studies | Tractability | Risk |
|---|---|---|---|---|
| **S — peer-reviewed supplement (Data C)** | Curated persulfidation site/protein tables from the paper | 006140, 024061, 039999 (+035795 fallback) | High (small files, public) | Supplement may be **protein-level only** → no `site_ms` |
| **M — registered mzid + method mapping** | Already-downloaded PXD035795 mzid + reviewed tag→site rule | 035795 | High (in hand) | Mass/tag→persulfidation mapping must survive review; CV term is "unknown modification" |
| **H — parse large search outputs** | Mascot `.dat`, MaxQuant `txt` zip member | 006140, 024061 | Low (GB-scale, parser-heavy) | Engineering cost; keep as cross-check only |

**Recommended A2 sequence:** attempt **M** for PXD035795 first (lowest new-download
cost, evidence already in hand), then **S** for one additional study
(PXD024061 preferred for its CC-BY license), to reach the Gate-1 minimum of two
distinct site-level studies. Treat **H** as validation, not the critical path.

> **Superseded by §5 empirical findings:** this sequence assumed supplements
> would be site-level. Inspection shows they are not; the practical outlook is
> now Phase Z. The route ranking is retained as the *method*, but its optimistic
> premise no longer holds.

---

## 5. Empirical inspection log (2026-07-21) — the actual answer to "is there site-level evidence?"

This section records a **real download-and-inspect pass** over every public
output, not a plan. All items were fetched over open channels and hashed.

### 5.1 What was acquired (open-access, hashed)

| Study | Source | Channel | SHA256 (bundle) | Bytes |
|---|---|---|---|---|
| PXD035795 | `peptides_1_1_0.mzid.gz` (in hand) | PRIDE FTP | `62105df…4d1b5fe` | 184,107 |
| PXD035795 | `Persulfidation_protects.pdf` (method) | IDUS Seville | `90b8c49…c43a554` | 11,657,462 |
| PXD006140 | PMC5853657 supplement bundle | Europe PMC | `3f48786…925d05` | 30,852,740 |
| PXD024061 | PMC8064375 supplement bundle | Europe PMC | `97fc78f…5b0066` | 4,110,151 |
| PXD039999 | PMC10433926 supplement bundle | Europe PMC | `d1eea46…5a1570bd` | 961,629 |

PMC map: PXD006140→PMC5853657, PXD024061→PMC8064375, PXD039999→PMC10433926
(via NCBI eLink from the registered PubMed IDs).

### 5.2 Per-study verdict (from real headers + method text)

| Study | Public **per-Cys site** table? | What the public output actually contains |
|---|---|---|
| **PXD006140** | **YES — via Dataset S3** (see §5.5) | The *summary* tables (S6 "S-sulfhydrated in WT/des1", S8 "Persulfidated proteins") are protein-level, BUT **Dataset S3** (`ident_peptides`) records the authors' **own custom persulfidation PTMs** — `MOD:99998 Sulfide`, `MOD:99997/99996 CN-Biotin-Sulfide` on **C** — with peptide-internal position (`MOD:99998 C14`) and protein coordinates (`inferredCoords`). This yields **residue-resolved, author-designated S-sulfhydration sites at FDR<1%**. This is genuine site-level evidence, not a reinterpretation. |
| **PXD024061** | **YES — via `txt_persulfproject.zip`** (see §5.6) | The MDPI *supplement* is protein-level only, BUT the deposited **MaxQuant search output** contains dedicated PTM-site tables `Sulfide(C)Sites.txt` (82 rows) and `CianoBiotin(C)Sites.txt` (6 rows) with `Positions within proteins`, `Amino acid=C`, `Localization prob`, `Reverse`/`Potential contaminant` flags — **~76 class-I (loc prob ≥0.75) residue-resolved persulfidation sites**. Same lab's tag-switch chemistry as PXD006140. |
| **PXD035795** | **NO** (honest) | mzid has **61 Cys-localized** modifications (NBF_C×44, DCP×14, dcp-ac×3) with peptide+location+protein. But Methods define **DCP (+168.08) = "sulfenic acid-reactive"** and **NBF (+163.00) = general Cys/Lys/Arg label**; the persulfidation-specific tag **Daz-2/Cy5 was used only for in-gel fluorescence, not LC–MS/MS**. So no MS tag is an author-sanctioned persulfidation-site mark. |
| **PXD039999** | **NO** | Dataset S1 (titled "Persulfida…") is protein-level: `Protein ID / #Peptides / PTM / [Control|Drought] intensities`, with **PTM column = "Carbamidomethylation"** (generic alkylation). Dataset S2 = differential protein abundance. No Cys-site column. |

### 5.3 Cross-study conclusion (revised after the §5.5 and §5.6 deep-dives)

**2 of 4 public studies yield genuine, author-designated, residue-resolved
persulfidation site evidence:**
- **PXD006140** — Dataset S3, ~356 sites (§5.5);
- **PXD024061** — MaxQuant `Sulfide(C)`/`CianoBiotin(C)` site tables in the
  deposited search output, ~76 class-I sites (§5.6).

The other two do not: PXD039999 publishes only a protein-level persulfidome;
PXD035795's MS tags (DCP/NBF) are sulfenic-acid/general rather than an
author-sanctioned persulfidation-site mark, and its New Phytologist article is not
in PMC / not OA.

**Implication for Gate 1 (GO is now plausible from public data):** the ≥2
distinct-site-level-study minimum **can be met** with PXD006140 + PXD024061. The
STOP was a coverage/selection matter (Blocker Types A/B), not a Type-D public-data
ceiling. Recommended path is now **forward to Phase B** (real PU benchmark) once
A2 lands both parsers under TDD — *not* Phase Z.

**Mandatory limitation to record in every downstream claim:** both confirmed
studies are from the **same laboratory (Romero/Gotor, Seville) using the same
tag-switch chemistry**, in the same species (Arabidopsis). Leave-study-out is
valid, but the two studies are **not chemically independent**; method- and
species-diversity are limited. A truly independent chemistry/lab/species remains a
priority **collaboration data request** to the Zhang lab (Codex §10) and keeps
`leave-species-out` (Split C) meaningful only once tomato/other data arrive.

### 5.4 Honest caveats on this pass

- The PXD024061/PXD039999 "protein-level only" verdicts are from sheet **headers**;
  I did not find a hidden peptide-site sheet in their bundles, but a full
  cell-level sweep was not performed.
- The PXD035795 New Phytologist supplement was **not** inspected in this pass — it
  is the top open lead for study #2 and must be fetched next.
- The PXD006140 site extraction (§5.5) is reconnaissance, not the production
  parser: it must be re-implemented under TDD in A2 with the exact conflict/decoy/
  isoform handling below.
- None of these bundles is yet registered in `supplementary_sources.tsv`; §6
  specifies that registration for whichever bundle A2 elects to parse.

### 5.5 Deep-dive — PXD006140 Dataset S3 **confirms** extractable site evidence

`erx294_suppl_supplementary_data_set_s3.xlsx` (sheet `ident_peptides`, 58,570 rows,
titled *"Quantitative comparison of S-sulfhydration patterns in wild type and des1
plants (FDR <1%)"*).

**Why it is author-designated persulfidation (not reinterpretation):** the
`Search Settings` sheet declares **custom variable modifications on Cys** built for
the tag-switch persulfidation assay:

| Code | Type | Site | Meaning | Role |
|---|---|---|---|---|
| `MOD:99998` | PTM | C | **Sulfide** (–SSH) | **persulfidation** |
| `MOD:99997` | PTM | C | **CN-Biotin-Sulfide** | **persulfidation (labeled)** |
| `MOD:99996` | PTM | C | CN-Biotin-Na-Sulfide | **persulfidation (labeled)** |
| `MOD:99999` | chemical | C | MSBT | free-thiol **control** (exclude) |
| `MOD:00110` | chemical | C | L-cysteine methyl disulfide | artifact (exclude) |

**Localization is explicit:** `mods` encodes `MOD:<code> <residue><pep_pos>`
(e.g. `MOD:99998 C14`) and `inferredCoords` gives protein span(s)
(e.g. `P27140-2{264-279}`), so protein position = span_start + pep_pos − 1.

**Yield (persulfidation mods 99996/99997/99998 on C, decoys dropped):**
- 435 non-decoy persulfidation peptide-spectra (305 decoy spectra correctly excluded);
- **356 distinct (protein, Cys-position) sites**, of which **323 from unambiguous
  single-protein peptides** and **51 spectra multi-mapped across isoforms** →
  belong in a conflict/isoform table, never silently resolved (matches Codex rules);
- occurrences: `MOD:99996`×454, `MOD:99998`×296, `MOD:99997`×47.

**Coordinate cross-check against real UniProt (reproducibility):**
- `ELGVPIVMHDYLTGGFTANTSLSHYCR` `MOD:99998 C26` `O03042{259-285}` → **O03042 Cys284 = C ✓**
- `GILAMDESNATCGK` `MOD:99998 C12` `Q944G9{68-81}` → **Q944G9 Cys79 = C ✓**

**Verdict:** PXD006140 Dataset S3 is a **Data-level-C, author-designated,
residue-resolved persulfidation site source (~356 sites, FDR<1%)**. It is the
first confirmed `site_*`-eligible study and the natural A2 Route-S/H target.

### 5.6 Deep-dive — PXD024061 **confirms** study #2 via MaxQuant site tables

The MDPI supplement is protein-level, but the **deposited MaxQuant search output**
`txt_persulfproject.zip` (registered file, PRIDE SHA1 `d0cbe145…`, 2.42 GB) was
inspected **without full download**: an HTTP Range read of the zip's central
directory listed 22 members, including dedicated MaxQuant **PTM-site tables** built
from the same tag-switch chemistry as PXD006140:

| Member | Rows | Meaning | Extracted-member SHA256 | Bytes |
|---|---|---|---|---|
| `Sulfide(C)Sites.txt` | 82 | persulfidation (–SSH) sites on C | `e9a1abfb…891003d8` | 151,666 |
| `CianoBiotin(C)Sites.txt` | 6 | CN-Biotin-Sulfide (labeled) sites on C | `a98a49e9…365d00fb` | 17,083 |
| `MSBT(C)Sites.txt` | — | free-thiol **control** (exclude) | — | 274,479 |

Only the two persulfidation members (169 KB total) were range-extracted and hashed;
the 2.42 GB archive was not materialized.

**Content (MaxQuant `…Sites.txt` schema, 363 columns):** `Proteins`,
`Leading proteins`, `Positions within proteins`, `Amino acid`, `Localization prob`,
per-sample localization probs, `Sequence window`, `Reverse`, `Potential contaminant`.

**Yield after decoy/contaminant removal + MaxQuant class-I filter (loc prob ≥0.75):**
- `Sulfide(C)Sites`: 82 → 80 non-decoy/contaminant → **70 class-I distinct
  (protein, Cys) sites** (all `Amino acid=C`); the 10 excluded are low-probability
  (e.g. adjacent C1019/C1021 in one CWC motif, prob 0.5) → lower-confidence tier,
  flagged not dropped silently.
- `CianoBiotin(C)Sites`: 6 → **6 class-I sites** (all C, prob ≥0.9997).
- Combined ≈ **76 class-I persulfidation sites**.

**Verdict:** PXD024061's deposited search output is a **Data-level-B/C,
author-designated, residue-resolved persulfidation site source (~76 class-I sites)**
— the second `site_*`-eligible study, satisfying Gate 1's two-study minimum. A2 must
register the two extracted members with their own SHA256 (parent archive SHA1 as
container provenance) and re-derive sites under TDD with the localization-prob
threshold, decoy/contaminant exclusion, and UniProt residue check.

---

## 6. Provenance registration plan

Reuse the existing hash-first, temp-file download policy and the
`evidence_methods.tsv` precedent. Add, per newly acquired supplementary file, a row
to a new registry table **`data/registry/supplementary_sources.tsv`** with columns:

```
study_accession
publication_doi
supplement_id          # e.g. "Table S3" / official supplementary filename
repository             # publisher | DOI resolver | institutional repo (e.g. IDUS)
official_url
local_path             # ../raw/supplements/<ACC>/<file>
retrieved_at
size_bytes
sha256
data_level             # C (peer-reviewed supplement) by default
license_or_usage       # observed at download; do NOT assume CC-BY for non-MDPI
scientific_use         # e.g. site_evidence_candidate | protein_level_context
parser_version
```

Rules (inherited from `AGENTS.md` / Codex §6):
- Download only after registering the expected source; verify SHA256 before use;
  a failed download must fail closed — never substitute a placeholder file.
- Manually review each supplement once and record the audit outcome (per Codex
  §3.3 step 4) before any site is emitted.
- Multi-GB Route-H members: register the **extracted member's** own SHA256, and
  record the parent archive's SHA256 as the container provenance.

---

## 7. Decisive open questions → Gate-1 decision matrix

The whole GO/STOP fork hinges on one empirical question per study, which **A2 must
resolve by inspecting the real file — not assume**:

| Study | Open question A2 must answer | If site-level | If protein-level only |
|---|---|---|---|
| PXD035795 | Does a registered method mapping justify a mass/tag on C as persulfidation? | → `site_ms` records | stays `unresolved` |
| PXD024061 | Does the CC-BY supplement list specific persulfidated Cys positions? | → `site_*` (Data C) | protein-level context only |
| PXD006140 | Does the 2017 supplement give Cys positions (vs proteins)? | → `site_*` (Data C) | protein-level context only |
| PXD039999 | Does the 2023 supplement give Cys positions? | → `site_*` (Data C) | protein-level context only |

**Gate-1 resolution (evaluated in A3):**
- **GO** if ≥ 2 distinct studies yield verified `site_*` evidence → proceed to
  Phase B (real PU benchmark).
- **STOP (persisting)** if fewer than 2 do → execute **Phase Z**: reframe the
  deliverable as *"植物硫巯基化数据资源与系统性证据审计"* (data-resource + evidence
  audit), and use this map as the backbone of the "what we can and cannot obtain
  from public data" argument and the collaboration data request to the Zhang lab.

---

## 8. A1 acceptance criteria (this document)

- [ ] Every registered study has an enumerated candidate-source table grounded in
      `files.tsv` and the verified DOI.
- [ ] Each candidate is labeled with expected data level and a site-vs-protein
      open question.
- [ ] The three acquisition routes are ranked with an explicit A2 sequence.
- [ ] A provenance schema for new supplementary sources is specified, reusing the
      `evidence_methods.tsv` precedent and the hash-first download policy.
- [ ] The Gate-1 decision matrix ties each study's open question to GO vs Phase Z.
- [ ] No parser, label, or evidence-class change is introduced by this cycle.

## 9. Handoff to A2 (first executable TDD cycle) — revised after §5.5

**Recommended first A2 target (now the strongest real evidence): PXD006140
Dataset S3 site parser.**

> Register `erx294_suppl_supplementary_data_set_s3.xlsx` (PMC5853657, Europe PMC)
> in `supplementary_sources.tsv` with its SHA256 and data-level C. Write a failing
> test, driven by a real micro-fixture cut from `ident_peptides`, that:
> - extracts only Cys carrying an author-designated persulfidation PTM
>   (`MOD:99998`/`99997`/`99996`), **excluding** `MOD:99999` MSBT and `MOD:00110`;
> - drops `is_decoy` rows;
> - computes protein position = `inferredCoords` span_start + pep_pos − 1 and
>   **verifies the residue is C against the registered UniProt sequence**
>   (reuse `test_site_coordinates_match_sequence`);
> - routes isoform multi-mapped spectra to a conflict table, never silently
>   picking the first;
> - emits `evidence_level=site_ms` only for validated, single-or-conflict-logged
>   sites; biological values unchanged, `source_manifest.json` recorded.
> Expected yield ≈ 356 sites (323 unambiguous). This makes PXD006140 the first
> `site_*`-eligible study.

**Second A2 target (satisfies Gate-1 ≥2 studies): PXD024061 MaxQuant site parser.**

> Register the two range-extracted members `Sulfide(C)Sites.txt`
> (SHA256 `e9a1abfb…`) and `CianoBiotin(C)Sites.txt` (SHA256 `a98a49e9…`) in
> `supplementary_sources.tsv` (parent `txt_persulfproject.zip` SHA1 `d0cbe145…` as
> container provenance; data-level B/C). Write a failing test, driven by a real
> micro-fixture, that:
> - keeps only `Amino acid=C` rows with `Reverse`≠`+` and `Potential contaminant`≠`+`;
> - applies the MaxQuant **class-I** threshold `Localization prob ≥ 0.75`, routing
>   lower-probability rows to a flagged lower-confidence tier (never dropped silently);
> - reads `Leading proteins` + `Positions within proteins` and **verifies the residue
>   is C against the registered UniProt sequence**;
> - emits `evidence_level=site_ms` for validated class-I sites.
> Expected yield ≈ 76 class-I sites. This makes PXD024061 the second
> `site_*`-eligible study → Gate 1 can flip to GO in A3.

*(The PXD035795 New Phytologist supplement is no longer on the critical path — the
article is not in PMC / not OA and study #2 is already secured. It remains an
optional third-study lead requiring institutional access.)*

> The earlier "PXD035795 mzid via Route M" cycle (module
> `src/plantpersulf/evidence/site_evidence.py` + its RED/GREEN test) remains valid
> as the **fail-closed gating mechanism**, but §5.5 shows the mzid tags are not an
> author-sanctioned persulfidation-site mark, so it should not be wired to emit
> real `site_ms` for PXD035795 without a genuine reviewed mapping. Keep it as the
> gate; do not use it to manufacture positives.

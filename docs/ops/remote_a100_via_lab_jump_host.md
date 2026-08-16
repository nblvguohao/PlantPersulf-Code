# 远程 GPU 服务器接入：本机 → 实验室主机 → A100（可跨项目复用）

> 首次跑通：2026-07-24，PlantPersulf-Code 项目，用于在无本地 Linux/WSL 环境下跑
> MMseqs2 真聚类（`protein_clusters_v2.tsv`，见 `phase_z_evidence_audit.md` §3.3）。
> 本文档与配套脚本 `scripts/remote/lab_a100_hop.py` 是**通用工具**，与
> PlantPersulf 项目本身无强绑定——复制这一个文件 + 两个凭据文件到任何其他项目
> 目录即可直接用。

## 1. 为什么需要两跳

A100 服务器（`100.112.165.109`）只信任来自实验室内网/Tailscale 已知节点的连接；
从任意外部机器直接 SSH 上去，即使密码正确，也会被拒绝（已实测验证：TCP 握手和
SSH banner 交换都正常，但 `Permission denied (publickey,password)`）。

实验室有一台 Windows 主机（`100.66.9.81`）常年在该内网里，且**已经配好了对 A100
的免密钥登录**（`~/.ssh/id_rsa` 已被 A100 授权）。所以唯一能走通的路径是：

```
本机 ──(SSH, 密码登录)──> 实验室 Windows 主机 ──(SSH, 已配好的密钥, 免密码)──> A100
```

两跳性质不同，都要注意：

- 实验室主机跑的是 **Windows OpenSSH Server**，远程执行环境是 `cmd.exe`，不是
  bash——发给它的命令必须是 cmd 语法。
- 转发给 A100 的命令是 **Linux bash 命令**，要被塞进一层 cmd.exe 转义里再发出去。
  两种转义语法嵌套极易出错（引号、反斜杠、管道符互相打架），所以本工具**把每条
  转发给 A100 的命令都 base64 编码**，彻底跳过转义问题。

## 2. 凭据文件（不进 git）

仓库根目录下两个纯文本文件，格式一致，均已在 `.gitignore` 里屏蔽：

**`LAB_HOST`**（跳板机）：
```
administrator@100.66.9.81
<密码>
G:\cc
```

**`A100`**（目标 GPU 服务器）：
```
user@100.112.165.109
<密码，实际未使用——目标机走密钥登录>
data/lgh
```

第三行是"这台机器上约定的工作目录"的备注，纯提示性，脚本不强制使用。

**换项目/换服务器复用时**：把 `scripts/remote/lab_a100_hop.py` 和这两个凭据文件
模板复制到新项目根目录，改内容即可，脚本本身不用改。

## 3. 工具：`scripts/remote/lab_a100_hop.py`

```bash
# 体检：两跳是否都通，GPU/磁盘/内存概况
python scripts/remote/lab_a100_hop.py check

# 在 A100 上跑一条 bash 命令
python scripts/remote/lab_a100_hop.py run "nvidia-smi -L"

# 本地文件 → A100（经跳板机中转，自动 SHA256 校验）
python scripts/remote/lab_a100_hop.py upload local_file.fasta /data/lgh/myproject/local_file.fasta

# A100 → 本地文件（同样自动校验）
python scripts/remote/lab_a100_hop.py download /data/lgh/myproject/result.tsv result.tsv
```

也可以当库用：

```python
from scripts.remote.lab_a100_hop import Hop
hop = Hop()
hop.run("mkdir -p /data/lgh/myproject")
hop.upload("data.fasta", "/data/lgh/myproject/data.fasta")
hop.run("cd /data/lgh/myproject && some_tool data.fasta")
hop.download("/data/lgh/myproject/out.tsv", "out.tsv")
hop.close()
```

依赖：仅 `paramiko`（`pip install paramiko`）。不依赖 Windows 自带的 `ssh.exe`
在本机上可用——所有 SSH 都是从 Python 里发起，第二跳则是**在跳板机上**调用它
自己的 `ssh`/`scp`。

## 4. 已知坑

### 4.1 Git Bash（MSYS）会偷偷改写看起来像绝对路径的参数

如果从 Git Bash 里跑这个脚本，任何形如 `/data/lgh/...` 的命令行参数会被 MSYS
自动"翻译"成 Windows 路径（比如变成 `C:/Program Files/Git/data/lgh/...`），
因为 MSYS 把它当成了本机的 POSIX 路径。这只发生在**远程路径参数**上——本地
路径参数其实需要这个转换（否则 Python 在 Windows 上收到的 `/tmp/x` 会被解析成
当前盘符下的 `\tmp\x`，同样是错的）。

**修复**：远程路径参数前面多写一个斜杠（`//data/lgh/...`），MSYS 看到 `//`
开头会当成 UNC 路径直接跳过转换：

```bash
# 本地路径：单斜杠，让 MSYS 正常转换
# 远程路径：双斜杠，阻止 MSYS 转换
python scripts/remote/lab_a100_hop.py upload /tmp/x.fasta //data/lgh/myproject/x.fasta
```

从 PowerShell / cmd.exe 跑則没有这个问题（没有 MSYS 这层），远程路径直接写单
斜杠即可。

### 4.2 Windows 主机的命令输出是 GBK 编码

跳板机是中文 Windows，`cmd.exe` 的默认输出编码是 GBK，不是 UTF-8。脚本内部已经
用 `decode('gbk', errors='replace')` 处理，直接用工具不会踩到这个坑；但如果自己
用 paramiko 手写连接，`stdout.read().decode('utf-8')` 遇到中文路径/日期会抛
`UnicodeDecodeError`，一定要用 `'gbk'`。

### 4.3 A100 上没有 sudo，没有 conda

目标服务器是多课题组共享机器，当前账号无密码 sudo（`sudo -n true` 会失败），
也没预装 conda/mamba。装工具走**用户空间静态二进制**（比如 MMseqs2 官方
`mmseqs-linux-avx2.tar.gz`），不要假设能装系统包。

### 4.4 `/data/lgh` 是多项目共享目录

服务器上这个目录下已经有一堆别的课题的文件夹（`GP_WAUTER_XGBoost`、
`plantomics-*`、`GWAS` 等）。**任何新工作都建一个专属子目录**
（如 `/data/lgh/<project_name>_<purpose>/`），不要把文件直接扔进 `/data/lgh/` 根。

### 4.5 跳板机上的中转文件要清理

`upload`/`download` 会在跳板机 `C:\Users\<lab_user>\_hop_transfer\` 下落一份
临时文件，用完自动删（脚本已处理）；如果脚本异常中断，记得手动检查这个目录，
不要在别人常用的机器上留垃圾。

## 5. 实测记录（首次跑通，供参考）

- 跳板机：`100.66.9.81`，Windows，`OpenSSH_for_Windows_8.1p1`。
- A100：`100.112.165.109`，Ubuntu 22.04（内核 6.8），双卡
  `NVIDIA A100-SXM4-80GB`，128 核，251GB 内存，`/data` 分区 15T（用 6.8T）。
- 用途：上传拟南芥参考蛋白组（54,646 条序列，28MB FASTA），跑
  `mmseqs easy-cluster --min-seq-id 0.3 -c 0.5 --cov-mode 0 --threads 32`，
  **6.9 秒**跑完（54,646 序列 → 13,967 簇），下载结果，SHA256 全程校验。
  详见 `docs/phase_z_evidence_audit.md` §3.3 / §5（MMseqs2 真聚类记录）。

## 6. 换到别的项目时要做的事

1. 复制 `scripts/remote/lab_a100_hop.py` 到新项目（原样，不用改）。
2. 在新项目根目录建 `LAB_HOST` 和 `<TARGET_SERVER_NAME>` 两个凭据文件
   （3 行格式：`user@host` / 密码 / 工作目录备注）。
3. 在新项目 `.gitignore` 里加上这两个文件名，避免误提交。
4. 如果目标服务器不是 A100 而是别的机器，把脚本里 `DEFAULT_TARGET_FILE` 指向
   的文件名换成新服务器的凭据文件名（或者调用时用
   `Hop(target_file=Path("OTHER_SERVER"))` 显式指定）。
5. `python scripts/remote/lab_a100_hop.py check` 先体检一遍，确认两跳都通。

"""Reusable two-hop remote execution: this machine -> lab Windows host -> A100 GPU server.

Why this exists
----------------
The A100 server only accepts SSH connections from the lab network (its sshd
rejects password auth for hosts outside it, even with the right password —
confirmed by testing directly from an off-network machine). The lab Windows
host sits inside that network and already has a passwordless SSH key
authorized on the A100 server. So the only working path from an arbitrary
outside machine is:

    this machine --(SSH, password auth)--> lab Windows host
                        --(SSH, key auth, already set up)--> A100 server

Both hops matter:
- The lab host runs OpenSSH Server for Windows (cmd.exe as the shell), so
  commands sent to it must be cmd-safe.
- The command *forwarded* to the A100 server is a Linux bash command, so it
  must survive being embedded inside a cmd.exe-quoted string. Nesting two
  incompatible quoting dialects is exactly the kind of thing that breaks
  silently, so every remote command run on A100 is base64-encoded before
  being sent — no quoting to get wrong, no encoding surprises.

Credential files (untracked, git-ignored, 3 lines each)
---------------------------------------------------------
    line 1: user@host
    line 2: password
    line 3: default working directory on that host (informational only)

Two files are expected at the repo root: ``LAB_HOST`` (the jump host) and
``A100`` (the target GPU server) — see ``.gitignore`` for why they must never
be committed. Copy this script and both credential files into any other
project that needs the same server; nothing here is PlantPersulf-specific.

CLI usage
---------
    python lab_a100_hop.py check
        Sanity-check both hops (whoami, nvidia-smi, disk/mem) and print a report.

    python lab_a100_hop.py run "<bash command>"
        Run a command on the A100 server, print stdout/stderr, return its exit code.

    python lab_a100_hop.py upload <local_path> <remote_path_on_a100>
        Copy a local file up to the A100 server (via the lab host), then
        SHA256-verify it landed intact.

    python lab_a100_hop.py download <remote_path_on_a100> <local_path>
        Copy a file back down from the A100 server, then SHA256-verify it.

Library usage
-------------
    from lab_a100_hop import Hop
    hop = Hop()                      # reads LAB_HOST / A100 from repo root
    hop.run("nvidia-smi -L")
    hop.upload("local.fasta", "/data/lgh/myproject/local.fasta")
    hop.download("/data/lgh/myproject/out.tsv", "out.tsv")
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import posixpath
import sys
import time
from dataclasses import dataclass
from pathlib import Path

try:
    import paramiko
except ImportError as exc:  # pragma: no cover
    raise SystemExit("paramiko is required: python -m pip install paramiko") from exc

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LAB_HOST_FILE = REPO_ROOT / "LAB_HOST"
DEFAULT_TARGET_FILE = REPO_ROOT / "A100"


@dataclass
class _Credential:
    username: str
    host: str
    password: str
    workdir: str


def _load_credential(path: Path) -> _Credential:
    if not path.exists():
        raise SystemExit(
            f"credential file not found: {path}\n"
            "expected 3 lines: user@host / password / default working dir"
        )
    lines = path.read_text(encoding="utf-8").splitlines()
    lines = [ln.strip() for ln in lines if ln.strip()]
    if len(lines) < 2:
        raise SystemExit(f"credential file {path} needs at least user@host + password")
    user_host, password = lines[0], lines[1]
    workdir = lines[2] if len(lines) > 2 else ""
    if "@" not in user_host:
        raise SystemExit(
            f"credential file {path} line 1 must be user@host, got: {user_host!r}"
        )
    username, host = user_host.split("@", 1)
    return _Credential(username=username, host=host, password=password, workdir=workdir)


class Hop:
    """Two-hop SSH session: local -> lab Windows jump host -> target Linux server.

    The jump host is reached with password auth (paramiko). The target server
    is reached by shelling out to the jump host's own `ssh`/`scp` client,
    which already carries a key authorized on the target — we never need the
    target's password.
    """

    def __init__(
        self,
        lab_host_file: Path = DEFAULT_LAB_HOST_FILE,
        target_file: Path = DEFAULT_TARGET_FILE,
        timeout: int = 20,
    ) -> None:
        self.lab = _load_credential(lab_host_file)
        self.target = _load_credential(target_file)
        self._client: paramiko.SSHClient | None = None
        self._timeout = timeout

    def _jump(self) -> paramiko.SSHClient:
        if self._client is None:
            client = paramiko.SSHClient()
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            client.connect(
                self.lab.host,
                username=self.lab.username,
                password=self.lab.password,
                look_for_keys=False,
                allow_agent=False,
                timeout=self._timeout,
            )
            self._client = client
        return self._client

    def _exec_on_lab_host(self, cmd: str, timeout: int = 60) -> tuple[int, str, str]:
        """Run a cmd.exe command on the lab host itself (not forwarded further)."""
        client = self._jump()
        stdin, stdout, stderr = client.exec_command(cmd, timeout=timeout)
        out = stdout.read().decode("gbk", errors="replace")
        err = stderr.read().decode("gbk", errors="replace")
        rc = stdout.channel.recv_exit_status()
        return rc, out, err

    def run(
        self, remote_cmd: str, timeout: int = 300, quiet: bool = False
    ) -> tuple[int, str, str]:
        """Run a bash command on the A100 target, base64-piped through the lab-host hop."""
        b64 = base64.b64encode(remote_cmd.encode("utf-8")).decode("ascii")
        hop_cmd = (
            f"ssh -o StrictHostKeyChecking=no {self.target.username}@{self.target.host} "
            f'"echo {b64} | base64 -d | bash"'
        )
        rc, out, err = self._exec_on_lab_host(hop_cmd, timeout=timeout)
        if not quiet:
            print(f">>> {remote_cmd}")
            if out:
                print(out)
            if err.strip():
                print("ERR:", err)
            print(f"[exit {rc}]\n")
        return rc, out, err

    def upload(
        self, local_path: str | Path, remote_path: str, verify: bool = True
    ) -> None:
        """Copy local_path -> A100:remote_path, staging through the lab host."""
        local_path = Path(local_path)
        if not local_path.exists():
            raise FileNotFoundError(local_path)

        local_sha = hashlib.sha256(local_path.read_bytes()).hexdigest()
        stage_dir = rf"C:\Users\{self.lab.username}\_hop_transfer"
        stage_path = f"{stage_dir}\\{local_path.name}"

        self._exec_on_lab_host(f'mkdir "{stage_dir}" 2>nul & echo ok')
        client = self._jump()
        sftp = client.open_sftp()
        t0 = time.time()
        sftp.put(str(local_path), stage_path)
        sftp.close()
        print(
            f"[upload] local -> lab host: {local_path.name} in {time.time() - t0:.1f}s"
        )

        remote_dir = posixpath.dirname(remote_path)
        if remote_dir:
            self.run(f"mkdir -p {remote_dir}", quiet=True)
        rc, out, err = self._exec_on_lab_host(
            f'scp -o StrictHostKeyChecking=no "{stage_path}" '
            f"{self.target.username}@{self.target.host}:{remote_path}",
            timeout=300,
        )
        self._exec_on_lab_host(f'del /f /q "{stage_path}" & echo cleaned')
        if rc != 0:
            raise RuntimeError(f"scp to target failed: {err or out}")

        if verify:
            rc, out, _ = self.run(f"sha256sum {remote_path}", quiet=True)
            remote_sha = out.split()[0] if out.strip() else ""
            if remote_sha != local_sha:
                raise RuntimeError(
                    f"SHA256 mismatch after upload: local={local_sha} remote={remote_sha}"
                )
            print(f"[upload] verified SHA256 {local_sha}")

    def download(
        self, remote_path: str, local_path: str | Path, verify: bool = True
    ) -> None:
        """Copy A100:remote_path -> local_path, staging through the lab host."""
        local_path = Path(local_path)
        local_path.parent.mkdir(parents=True, exist_ok=True)

        remote_sha = ""
        if verify:
            rc, out, _ = self.run(f"sha256sum {remote_path}", quiet=True)
            if rc == 0 and out.strip():
                remote_sha = out.split()[0]

        stage_dir = rf"C:\Users\{self.lab.username}\_hop_transfer"
        fname = posixpath.basename(remote_path)
        stage_path = f"{stage_dir}\\{fname}"
        self._exec_on_lab_host(f'mkdir "{stage_dir}" 2>nul & echo ok')
        rc, out, err = self._exec_on_lab_host(
            f"scp -o StrictHostKeyChecking=no "
            f'{self.target.username}@{self.target.host}:{remote_path} "{stage_path}"',
            timeout=300,
        )
        if rc != 0:
            raise RuntimeError(f"scp from target failed: {err or out}")

        client = self._jump()
        sftp = client.open_sftp()
        t0 = time.time()
        sftp.get(stage_path, str(local_path))
        sftp.close()
        print(
            f"[download] lab host -> local: {local_path.name} in {time.time() - t0:.1f}s"
        )
        self._exec_on_lab_host(f'del /f /q "{stage_path}" & echo cleaned')

        if verify and remote_sha:
            local_sha = hashlib.sha256(local_path.read_bytes()).hexdigest()
            if local_sha != remote_sha:
                raise RuntimeError(
                    f"SHA256 mismatch after download: remote={remote_sha} local={local_sha}"
                )
            print(f"[download] verified SHA256 {local_sha}")

    def close(self) -> None:
        if self._client is not None:
            self._client.close()
            self._client = None


def _cmd_check(hop: Hop) -> None:
    print(f"jump host : {hop.lab.username}@{hop.lab.host}")
    print(
        f"target    : {hop.target.username}@{hop.target.host}  (workdir hint: {hop.target.workdir})"
    )
    rc, out, _ = hop._exec_on_lab_host("whoami & hostname")
    print("[lab host]\n" + out)
    hop.run("whoami && hostname && uname -a")
    hop.run("nvidia-smi -L 2>&1 || echo 'no GPU visible'")
    hop.run("df -h . 2>&1")
    hop.run("free -h 2>&1")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("check", help="Sanity-check both hops")

    p_run = sub.add_parser("run", help="Run a bash command on the A100 target")
    p_run.add_argument("remote_cmd")
    p_run.add_argument("--timeout", type=int, default=300)

    p_up = sub.add_parser("upload", help="Copy a local file to the A100 target")
    p_up.add_argument("local_path")
    p_up.add_argument("remote_path")

    p_down = sub.add_parser("download", help="Copy a file from the A100 target")
    p_down.add_argument("remote_path")
    p_down.add_argument("local_path")

    args = parser.parse_args()
    hop = Hop()
    try:
        if args.command == "check":
            _cmd_check(hop)
            return 0
        if args.command == "run":
            rc, _, _ = hop.run(args.remote_cmd, timeout=args.timeout)
            return rc
        if args.command == "upload":
            hop.upload(args.local_path, args.remote_path)
            return 0
        if args.command == "download":
            hop.download(args.remote_path, args.local_path)
            return 0
    finally:
        hop.close()
    return 1


if __name__ == "__main__":
    sys.exit(main())

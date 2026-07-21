# Data reproduction

Raw and intermediate datasets are intentionally excluded from Git. The
machine-readable download list is
`data/registry/reproduction_downloads.tsv`. It contains the accession, official
HTTPS URL, raw-data-relative destination, expected byte counts, repository
checksum where available, SHA256, and the restricted scientific use assigned in
this project.

From the repository root, the following PowerShell snippet downloads every
registered input into `data/raw` and rejects any file whose SHA256 does not
match the frozen manifest:

```powershell
$rows = Import-Csv data/registry/reproduction_downloads.tsv -Delimiter "`t"
foreach ($row in $rows) {
    $destination = Join-Path data/raw $row.relative_path
    New-Item -ItemType Directory -Force -Path (Split-Path $destination) |
        Out-Null
    Invoke-WebRequest -Uri $row.official_url -OutFile $destination
    $actual = (Get-FileHash -Algorithm SHA256 $destination).Hash.ToLower()
    if ($actual -ne $row.sha256) {
        throw "SHA256 mismatch: $destination"
    }
}
```

The list is a provenance and acquisition aid, not a label definition. In
particular, nondetection is not converted into an experimental negative, and
no benchmark label, split, model result, or performance metric is introduced
by downloading these files.

$ErrorActionPreference = "Stop"
python (Join-Path $PSScriptRoot "run_demo.py") @args
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

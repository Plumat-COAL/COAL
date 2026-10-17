<#
.SYNOPSIS
  Build the COAL desktop app for Windows: venv -> tests -> PyInstaller one-folder build ->
  frozen-build checks -> Inno Setup installer (if Inno Setup 6 is installed).

.DESCRIPTION
  Run from anywhere in a PowerShell window (Windows 10/11, x64):

      powershell -ExecutionPolicy Bypass -File desktop\build.ps1

  Needs: Python 3.12 from python.org (the "py" launcher), internet access for pip the first
  time, and optionally Inno Setup 6.3+ (https://jrsoftware.org/isdl.php) for the installer.

  Results:
      desktop\dist\COAL\COAL.exe            portable one-folder app (zip this folder to share)
      desktop\dist\COAL-Setup-<ver>.exe     per-user installer (only if Inno Setup was found)

  Every step has a timeout; a step that hangs is stopped and the build fails with a message.
  Nothing is pushed or uploaded anywhere.

.PARAMETER Python
  Python launcher arguments used to create the venv (default: "-3.12" for "py -3.12").
.PARAMETER Clean
  Delete the venv, build and dist folders first.
.PARAMETER SkipTests
  Skip the unit tests before building.
.PARAMETER NoInstaller
  Don't run Inno Setup even if it is installed.
#>
param(
    [string]$Python = "-3.12",
    [switch]$Clean,
    [switch]$SkipTests,
    [switch]$NoInstaller
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$Desktop = $PSScriptRoot
$Repo = Split-Path -Parent $Desktop
$Venv = Join-Path $Desktop ".venv"
$Py = Join-Path $Venv "Scripts\python.exe"
$Dist = Join-Path $Desktop "dist"
$Work = Join-Path $Desktop "build"
$Exe = Join-Path $Dist "COAL\COAL.exe"

function Quote-Arg([string]$a) {
    # Windows command-line quoting for one argument (handles spaces and embedded quotes).
    if ($a -eq "") { return '""' }
    if ($a -notmatch '[\s"]') { return $a }
    $escaped = $a -replace '(\\*)"', '$1$1\"'
    $escaped = $escaped -replace '(\\+)$', '$1$1'
    return '"' + $escaped + '"'
}

function Invoke-Step {
    # Run a program with a timeout; throw if it fails or takes too long.
    param([string]$Title, [string]$File, [string[]]$Arguments, [int]$TimeoutSec = 600,
          [string]$WorkingDirectory = $Repo)
    Write-Host ""
    Write-Host "==> $Title" -ForegroundColor Cyan
    $argLine = ($Arguments | ForEach-Object { Quote-Arg $_ }) -join " "
    Write-Host "    $File $argLine"
    $p = Start-Process -FilePath $File -ArgumentList $argLine -WorkingDirectory $WorkingDirectory `
        -NoNewWindow -PassThru
    $null = $p.Handle  # keep the handle so ExitCode is available after exit
    if (-not $p.WaitForExit($TimeoutSec * 1000)) {
        Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
        throw "$Title did not finish within $TimeoutSec seconds and was stopped."
    }
    $p.WaitForExit()
    if ($p.ExitCode -ne 0) { throw "$Title failed (exit code $($p.ExitCode))." }
}

Push-Location $Repo
try {
    if ($Clean) {
        foreach ($d in @($Venv, $Dist, $Work)) {
            if (Test-Path $d) { Write-Host "Removing $d"; Remove-Item -Recurse -Force $d }
        }
    }

    # 1. Local venv (never the system Python), with pinned requirements.
    if (-not (Test-Path $Py)) {
        Invoke-Step "Create venv desktop\.venv" "py" @($Python, "-m", "venv", $Venv) 300
    }
    Invoke-Step "Upgrade pip" $Py @("-m", "pip", "install", "--upgrade", "pip") 300
    Invoke-Step "Install COAL + GUI + desktop requirements" $Py `
        @("-m", "pip", "install", "-e", ".[gui]", "-r", "desktop\requirements.txt") 1200

    # 2. Tests (engine, app, wrapper, smoke). The smoke test starts the real local app.
    if (-not $SkipTests) {
        Invoke-Step "Engine tests" $Py @("-m", "unittest", "discover", "-s", "tests") 600
        Invoke-Step "App tests" $Py @("-m", "unittest", "discover", "-s", "tests_app") 900
        Invoke-Step "Wrapper tests" $Py @("-m", "unittest", "desktop\test_wrapper.py") 300
        Invoke-Step "Smoke test (source)" $Py @("-m", "unittest", "desktop\test_smoke.py", "-v") 900
    }

    # 3. PyInstaller one-folder build.
    $version = (Select-String -Path "src\coal\__init__.py" -Pattern '__version__\s*=\s*"([^"]+)"').Matches[0].Groups[1].Value
    Write-Host "COAL version $version"
    Invoke-Step "PyInstaller build" $Py @("-m", "PyInstaller", "--noconfirm", "--clean",
        "desktop\coal_desktop.spec", "--distpath", $Dist, "--workpath", $Work) 1800
    if (-not (Test-Path $Exe)) { throw "PyInstaller finished but $Exe is missing." }

    # 4. Check the frozen build: engine dispatch, then the smoke test against COAL.exe.
    Invoke-Step "Frozen engine: COAL.exe -m coal --version" $Exe @("-m", "coal", "--version") 120
    if (-not $SkipTests) {
        $env:COAL_SMOKE_EXE = $Exe
        try {
            Invoke-Step "Smoke test (frozen COAL.exe)" $Py @("-m", "unittest", "desktop\test_smoke.py", "-v") 900
        } finally {
            Remove-Item Env:\COAL_SMOKE_EXE -ErrorAction SilentlyContinue
        }
    }

    # 5. Installer, if Inno Setup 6 is installed (per-user or machine-wide).
    $iscc = @(
        (Join-Path $env:LOCALAPPDATA "Programs\Inno Setup 6\ISCC.exe"),
        (Join-Path ${env:ProgramFiles(x86)} "Inno Setup 6\ISCC.exe"),
        (Join-Path $env:ProgramFiles "Inno Setup 6\ISCC.exe")
    ) | Where-Object { $_ -and (Test-Path $_) } | Select-Object -First 1
    if ($NoInstaller) {
        Write-Host "Skipping the installer (-NoInstaller)."
    } elseif ($iscc) {
        Invoke-Step "Inno Setup installer" $iscc @("/Qp", "/DAppVersion=$version", "desktop\installer.iss") 900
        Write-Host ""
        Write-Host "Installer: $(Join-Path $Dist "COAL-Setup-$version.exe")" -ForegroundColor Green
    } else {
        Write-Host ""
        Write-Host "Inno Setup 6 not found, so no installer was made. Install it from" -ForegroundColor Yellow
        Write-Host "https://jrsoftware.org/isdl.php (per-user install is fine), then run:" -ForegroundColor Yellow
        Write-Host "  & `"`$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe`" /DAppVersion=$version desktop\installer.iss"
    }

    Write-Host ""
    Write-Host "Portable app: $Exe" -ForegroundColor Green
    Write-Host "Next: follow desktop\SMOKE_TEST.md on a clean Windows machine."
}
finally {
    Pop-Location
}

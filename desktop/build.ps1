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
        # CHANGED(WIN-W7): was Stop-Process (top process only; pip/PyInstaller/test children kept
        # running); now the whole tree. See overnight/windows-app/CHANGES.md
        & taskkill.exe /T /F /PID $p.Id | Out-Null
        throw "$Title did not finish within $TimeoutSec seconds and was stopped."
    }
    $p.WaitForExit()
    if ($p.ExitCode -ne 0) { throw "$Title failed (exit code $($p.ExitCode))." }
}

function Invoke-Capture {
    # Run a program with a timeout and return its stdout+stderr (redirected to files, because a
    # windowed exe has no console to print to). Throws on timeout or non-zero exit.
    param([string]$Title, [string]$File, [string[]]$Arguments, [int]$TimeoutSec = 120)
    Write-Host ""
    Write-Host "==> $Title" -ForegroundColor Cyan
    $argLine = ($Arguments | ForEach-Object { Quote-Arg $_ }) -join " "
    $out = [System.IO.Path]::GetTempFileName()
    $err = [System.IO.Path]::GetTempFileName()
    try {
        $p = Start-Process -FilePath $File -ArgumentList $argLine -WorkingDirectory $Repo `
            -NoNewWindow -PassThru -RedirectStandardOutput $out -RedirectStandardError $err
        $null = $p.Handle
        if (-not $p.WaitForExit($TimeoutSec * 1000)) {
            & taskkill.exe /T /F /PID $p.Id | Out-Null
            throw "$Title did not finish within $TimeoutSec seconds and was stopped."
        }
        $p.WaitForExit()
        $text = [string](Get-Content -Raw $out) + [string](Get-Content -Raw $err)
        if ($p.ExitCode -ne 0) { throw "$Title failed (exit code $($p.ExitCode)): $text" }
        return $text
    } finally {
        Remove-Item $out, $err -ErrorAction SilentlyContinue
    }
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
    # CHANGED(WIN-W7): was: exit code only (a windowed exe's output is invisible); now the output
    # is captured and must name the version. See overnight/windows-app/CHANGES.md
    $v = Invoke-Capture "Frozen engine: COAL.exe -m coal --version" $Exe @("-m", "coal", "--version") 120
    if ($v -notmatch [regex]::Escape($version)) {
        throw "COAL.exe -m coal --version printed '$($v.Trim())', expected version $version."
    }
    Write-Host "    $($v.Trim())"
    if (-not $SkipTests) {
        $env:COAL_SMOKE_EXE = $Exe
        try {
            Invoke-Step "Smoke test (frozen COAL.exe)" $Py @("-m", "unittest", "desktop\test_smoke.py", "-v") 900
            # CHANGED(SHIP-19): new. The bundled C++ engine (when built) must be found by the frozen app and give
            # the Python engine's bytes.
            $cppExe = Join-Path $Dist "COAL\_internal\coal-cpp.exe"
            if (Test-Path $cppExe) {
                $a = Join-Path $env:TEMP "coal_build_cpp.fasta"
                $b = Join-Path $env:TEMP "coal_build_py.fasta"
                $env:COAL_ENGINE = "cpp"
                try {
                    $null = Invoke-Capture "Frozen engine, C++ (COAL_ENGINE=cpp)" $Exe `
                        @("-m", "coal", "examples\clean\protein.aln.fasta", "examples\clean\dna.fasta", "-o", $a) 120
                } finally {
                    Remove-Item Env:\COAL_ENGINE -ErrorAction SilentlyContinue
                }
                Invoke-Step "Python engine (reference)" $Py `
                    @("-m", "coal", "examples\clean\protein.aln.fasta", "examples\clean\dna.fasta", "-o", $b) 120
                if ((Get-FileHash $a).Hash -ne (Get-FileHash $b).Hash) { throw "bundled coal-cpp output differs from Python" }
                Write-Host "    bundled coal-cpp: found by the app, output byte-identical to Python"
            } else {
                Write-Host "No coal-cpp.exe in the bundle (build it first: cmake -S . -B build; cmake --build build --config Release)." -ForegroundColor Yellow
            }
            # CHANGED(SHIP-18): new. What Claude Desktop starts after "Connect COAL to Claude Desktop": the frozen
            # app's stdio MCP server must answer and run a real alignment.
            Invoke-Step "MCP for Claude Desktop (frozen COAL.exe -m coal_app --mcp)" $Py `
                @("-m", "unittest", "-v", "tests_app.test_tools.TestFrozenAppMcp") 600
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

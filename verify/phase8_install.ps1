# Phase 8: uninstall the installed COAL 0.5.0, reinstall from the release file with the two optional tasks, check,
# uninstall again, then reinstall with the defaults and reconnect Claude Desktop (the state this machine started in).
# Silent installs skip the [Run] page (skipifsilent), so --connect-claude is run by hand where the wizard would.
$ErrorActionPreference = "Continue"
$out = "D:\coal-stress\verify\p8"
$setup = "D:\coal-handoff-kit\release\COAL-Setup-0.5.0.exe"
$app = "$env:LOCALAPPDATA\Programs\COAL"
$exe = "$app\COAL.exe"
$claude = "$env:APPDATA\Claude\claude_desktop_config.json"
$lnk = "$env:APPDATA\Microsoft\Windows\Start Menu\Programs\COAL.lnk"
$dlnk = [Environment]::GetFolderPath("Desktop") + "\COAL.lnk"
$checks = @()
function Check($name, $ok, $detail) { $script:checks += [pscustomobject]@{ name = $name; ok = [bool]$ok; detail = "$detail" }; if ($ok) { "OK   $name" } else { "FINDING-CANDIDATE $name : $detail" } }
function Uninst() { return Get-ItemProperty "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\*" -ErrorAction SilentlyContinue | Where-Object { $_.DisplayName -eq "COAL" } }
function OpenWith() { return @(".aln", ".fasta", ".fa", ".fas", ".faa", ".fna") | Where-Object { Test-Path "HKCU:\Software\Classes\$_\OpenWithProgids" -and (Get-Item "HKCU:\Software\Classes\$_\OpenWithProgids").Property -contains "COAL.SequenceFile" } }
function ProgId() { return Test-Path "HKCU:\Software\Classes\COAL.SequenceFile" }
function Snapshot($label) {
  $s = [ordered]@{ label = $label; app_dir = (Test-Path $app); exe = (Test-Path $exe); start_menu = (Test-Path $lnk); desktop = (Test-Path $dlnk)
    uninstall_entry = ($null -ne (Uninst)); openwith = @(OpenWith); progid = (ProgId); claude = $(if (Test-Path $claude) { Get-Content $claude -Raw } else { "ABSENT" })
    claude_bak = (Test-Path "$claude.bak-coal"); temp_gradio = (Test-Path "$env:TEMP\gradio"); temp_runs = (Test-Path "$env:TEMP\coal_app_runs"); data_dir = (Test-Path "$env:LOCALAPPDATA\COAL")
    fasta_default = (Get-ItemProperty "HKCU:\Software\Classes\.fasta" -ErrorAction SilentlyContinue).'(default)'
    fasta_userchoice = (Get-ItemProperty "HKCU:\Software\Microsoft\Windows\CurrentVersion\Explorer\FileExts\.fasta\UserChoice" -ErrorAction SilentlyContinue).ProgId }
  $s | ConvertTo-Json -Compress | Out-File "$out\state_$label.json" -Encoding utf8
  return $s
}
# the Claude config as the installer left it (coal only); add another server so uninstall can be checked for keeping it
$orig = Get-Content $claude -Raw
$before = Snapshot "0_before"
"--- state before: " + ($before | ConvertTo-Json -Compress)
$withOther = '{ "mcpServers": { "other": { "command": "C:\\x\\other.exe" }, "coal": { "command": "' + ($exe -replace '\\', '\\') + '", "args": ["-m", "coal_app", "--mcp"] } }, "theme": "dark" }'
Set-Content $claude $withOther -Encoding utf8
# 1. uninstall
$t = Measure-Command { $p = Start-Process "$app\unins000.exe" -ArgumentList "/VERYSILENT", "/NORESTART" -Wait -PassThru }
Check "uninstall (silent) exits 0" ($p.ExitCode -eq 0) "exit $($p.ExitCode), $([int]$t.TotalSeconds) s"
Start-Sleep -Seconds 2
$s1 = Snapshot "1_after_uninstall"
Check "uninstall removes the app folder" (-not $s1.app_dir) $s1.app_dir
Check "uninstall removes the Start Menu shortcut" (-not $s1.start_menu) $s1.start_menu
Check "uninstall removes the Add/Remove entry" (-not $s1.uninstall_entry) $s1.uninstall_entry
$c = $s1.claude | ConvertFrom-Json
Check "uninstall removes the coal entry from claude_desktop_config.json" ($null -eq $c.mcpServers.coal) $s1.claude
Check "uninstall keeps the other server and other keys" ($c.mcpServers.other.command -eq "C:\x\other.exe" -and $c.theme -eq "dark") $s1.claude
Check "uninstall leaves %TEMP%\gradio (documented)" $s1.temp_gradio $s1.temp_gradio
Check "uninstall leaves the per-user data folder %LOCALAPPDATA%\COAL (logs; not documented either way)" $true $s1.data_dir
Check "uninstall leaves .fasta's default program and UserChoice as they were" ($s1.fasta_default -eq $before.fasta_default -and $s1.fasta_userchoice -eq $before.fasta_userchoice) "$($s1.fasta_default) / $($s1.fasta_userchoice)"
# 2. reinstall with both optional tasks
$t = Measure-Command { $p = Start-Process $setup -ArgumentList "/VERYSILENT", "/NORESTART", "/TASKS=desktopicon,openwith", "/LOG=$out\install_tasks.log" -Wait -PassThru }
Check "install with /TASKS=desktopicon,openwith exits 0 (no UAC: per-user)" ($p.ExitCode -eq 0) "exit $($p.ExitCode), $([int]$t.TotalSeconds) s"
$s2 = Snapshot "2_after_install_tasks"
Check "install: app folder, COAL.exe, Start Menu shortcut, Add/Remove entry" ($s2.app_dir -and $s2.exe -and $s2.start_menu -and $s2.uninstall_entry) ($s2 | ConvertTo-Json -Compress)
Check "desktop shortcut task creates COAL.lnk on the Desktop" $s2.desktop $dlnk
Check "Open-with task registers COAL.SequenceFile for .aln .fasta .fa .fas .faa .fna" ($s2.progid -and $s2.openwith.Count -eq 6) ($s2.openwith -join ",")
Check "Open-with task leaves .fasta's default program and UserChoice unchanged" ($s2.fasta_default -eq $before.fasta_default -and $s2.fasta_userchoice -eq $before.fasta_userchoice) "$($s2.fasta_default) / $($s2.fasta_userchoice)"
$u = Uninst; Check "Add/Remove entry says COAL 0.5.0 in %LOCALAPPDATA%\Programs\COAL" ($u.DisplayVersion -eq "0.5.0" -and $u.InstallLocation -like "*\AppData\Local\Programs\COAL*") "$($u.DisplayVersion) $($u.InstallLocation)"
$log = Get-Content "$out\install_tasks.log" -Raw
Check "install log: no 'Error' / 'Failed' lines" (-not ($log -match "(?m)^.*(Error|Failed).*$")) (($log -split "`n" | Select-String "Error|Failed" | Select-Object -First 3) -join " | ")
$admin = ($log -match "IsAdmin: No" -or $log -match "Administrative install: No" -or $log -match "Non administrative install mode")
Check "install log says a non-administrative (per-user) install" $admin (($log -split "`n" | Select-String -Pattern "admin" -CaseSensitive:$false | Select-Object -First 3) -join " | ")
# the Open-with verb opens COAL with the file (check the command string only; the window itself is phase 8's Playwright part)
$cmd = (Get-ItemProperty "HKCU:\Software\Classes\COAL.SequenceFile\shell\open\command" -ErrorAction SilentlyContinue).'(default)'
Check "Open-with command is COAL.exe ""%1""" ($cmd -eq "`"$exe`" `"%1`"") $cmd
# the [Run] page is skipped by a silent install: run its first entry by hand against the config with another server (as the wizard would)
Set-Content $claude '{ "mcpServers": { "other": { "command": "C:\\x\\other.exe" } }, "theme": "dark" }' -Encoding utf8
Remove-Item "$claude.bak-coal" -ErrorAction SilentlyContinue
$p = Start-Process $exe -ArgumentList "--connect-claude", "--quiet" -Wait -PassThru
$c = Get-Content $claude -Raw | ConvertFrom-Json
Check "wizard's Connect step (COAL.exe --connect-claude): coal entry added, other server and keys kept, .bak-coal written" ($p.ExitCode -eq 0 -and $c.mcpServers.coal.command -eq $exe -and $c.mcpServers.other.command -eq "C:\x\other.exe" -and $c.theme -eq "dark" -and (Test-Path "$claude.bak-coal")) (Get-Content $claude -Raw)
# 3. uninstall again: the tasks' entries must go too
$p = Start-Process "$app\unins000.exe" -ArgumentList "/VERYSILENT", "/NORESTART" -Wait -PassThru
$s3 = Snapshot "3_after_uninstall_2"
Check "second uninstall removes the desktop shortcut, the Open-with entries and the ProgId" (-not $s3.desktop -and $s3.openwith.Count -eq 0 -and -not $s3.progid -and -not $s3.app_dir) ($s3 | ConvertTo-Json -Compress)
$c = $s3.claude | ConvertFrom-Json
Check "second uninstall removes coal from Claude's config, keeps the other server" ($null -eq $c.mcpServers.coal -and $c.mcpServers.other.command -eq "C:\x\other.exe") $s3.claude
Check "uninstall does not restore the .bak-coal copy over the current file (keeps the user's later edits)" ($c.theme -eq "dark") $s3.claude
# 4. final: defaults, reconnect Claude, restore the original config content (coal entry only, as the installer made it)
$p = Start-Process $setup -ArgumentList "/VERYSILENT", "/NORESTART", "/LOG=$out\install_final.log" -Wait -PassThru
Check "final install (defaults) exits 0" ($p.ExitCode -eq 0) $p.ExitCode
Remove-Item "$claude.bak-coal" -ErrorAction SilentlyContinue
Set-Content $claude $orig -Encoding utf8 -NoNewline
$s4 = Snapshot "4_final"
Check "final state equals the state before (app, shortcut, no desktop shortcut, no Open-with, Claude entry)" ($s4.app_dir -and $s4.start_menu -and -not $s4.desktop -and $s4.openwith.Count -eq 0 -and $s4.uninstall_entry -and ($s4.claude -replace '\s', '') -eq ($before.claude -replace '\s', '')) ($s4 | ConvertTo-Json -Compress)
$p = Start-Process $exe -ArgumentList "--version" -Wait -PassThru -NoNewWindow -RedirectStandardOutput "$out\final_version.txt"
Check "final COAL.exe --version" ((Get-Content "$out\final_version.txt" -Raw) -match "0\.5\.0") (Get-Content "$out\final_version.txt" -Raw)
$checks | ConvertTo-Json | Out-File "$out\install_results.json" -Encoding utf8
"$($checks.Count) checks, $(($checks | Where-Object { -not $_.ok }).Count) finding candidates"

; Inno Setup 6 script: per-user COAL installer (no administrator rights needed).
;
; Build (after desktop\build.ps1 has made desktop\dist\COAL\), from the repo root:
;   "%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe" /DAppVersion=0.3.0 desktop\installer.iss
; (or "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" for a machine-wide Inno Setup install).
; Output: desktop\dist\COAL-Setup-<version>.exe
;
; Needs Inno Setup 6.3 or newer (x64compatible). NOT BUILT YET: written on a Linux machine
; without Inno Setup (see overnight/windows-app/LOG.md).
;
; Decisions (DESIGN.md 5.4, 5.5): installs to %LOCALAPPDATA%\Programs\COAL; one Start menu
; shortcut; desktop shortcut and "Open with" entries are optional and off by default; the
; default program for .fasta/.aln is never changed; uninstall removes the program, shortcuts,
; registry entries and COAL's own app data, never results; no update check.

#ifndef AppVersion
  #define AppVersion "0.3.0"
#endif
#define AppName "COAL"
#define AppExe "COAL.exe"
#define ProgId "COAL.SequenceFile"

[Setup]
; Same AppId for every version, so a new installer upgrades in place. Never change it.
AppId={{E1184717-7056-4AF1-9288-5AB3BC967E81}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=COAL contributors
AppPublisherURL=https://github.com/chloe-newt/COAL
AppSupportURL=https://github.com/chloe-newt/COAL/issues
AppComments=Codon alignment from a protein alignment and coding DNA. Works offline.
; Per-user install: no UAC prompt. {autopf} is %LOCALAPPDATA%\Programs in this mode.
PrivilegesRequired=lowest
DefaultDirName={autopf}\{#AppName}
DisableDirPage=yes
DisableProgramGroupPage=yes
DisableReadyPage=yes
DisableWelcomePage=yes
; The wrapper holds this named mutex while running: setup and uninstall ask to close COAL first.
AppMutex=COAL-desktop-running
CloseApplications=yes
RestartApplications=no
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
WizardStyle=modern
SetupIconFile=coal.ico
UninstallDisplayIcon={app}\{#AppExe}
UninstallDisplayName={#AppName}
OutputDir=dist
OutputBaseFilename=COAL-Setup-{#AppVersion}
Compression=lzma2/max
SolidCompression=yes
; Tell Explorer to refresh "Open with" menus, only when that task was chosen.
; CHANGED(WIN-W8): was yes; now a boolean expression (Inno Setup 6 allows one here; not compiled
; yet: if ISCC rejects it, go back to yes). See overnight/windows-app/CHANGES.md
ChangesAssociations=WizardIsTaskSelected('openwith')
SetupLogging=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Messages]
SelectTasksLabel2=COAL installs for you only and needs no administrator rights. Choose any extras, then click Install.

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; Flags: unchecked
; CHANGED(WIN-W4): was: offered "Open with" with no caveat; now says what it does today. See overnight/windows-app/CHANGES.md
Name: "openwith"; Description: "Add COAL to ""Open with"" for .aln and .fasta files (opens COAL and tells you which box to drag the file onto; it can't load the file by itself yet; your default programs don't change)"; Flags: unchecked

[InstallDelete]
; CHANGED(WIN-W8): new. An upgrade first removes the old bundle, so files a newer build no longer
; ships don't linger and get imported. See overnight/windows-app/CHANGES.md
Type: filesandordirs; Name: "{app}\_internal"

[Files]
Source: "dist\COAL\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\LICENSE"; DestDir: "{app}"; DestName: "LICENSE.txt"; Flags: ignoreversion
Source: "README.md"; DestDir: "{app}"; DestName: "README.md"; Flags: ignoreversion skipifsourcedoesntexist

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\{#AppExe}"; Comment: "Codon alignment, offline"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Comment: "Codon alignment, offline"; Tasks: desktopicon

[Registry]
; "Open with" only (OpenWithProgids), per user (HKA = HKCU in a per-user install).
Root: HKA; Subkey: "Software\Classes\{#ProgId}"; ValueType: string; ValueName: ""; ValueData: "Sequence file"; Flags: uninsdeletekey; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\{#ProgId}\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\{#AppExe},0"; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\{#ProgId}\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#AppExe}"" ""%1"""; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\Applications\{#AppExe}"; ValueType: string; ValueName: "FriendlyAppName"; ValueData: "{#AppName}"; Flags: uninsdeletekey; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.aln\OpenWithProgids"; ValueType: string; ValueName: "{#ProgId}"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.fasta\OpenWithProgids"; ValueType: string; ValueName: "{#ProgId}"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.fa\OpenWithProgids"; ValueType: string; ValueName: "{#ProgId}"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.fas\OpenWithProgids"; ValueType: string; ValueName: "{#ProgId}"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.faa\OpenWithProgids"; ValueType: string; ValueName: "{#ProgId}"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.fna\OpenWithProgids"; ValueType: string; ValueName: "{#ProgId}"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith

[Run]
Filename: "{app}\{#AppExe}"; Description: "Open COAL now"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; COAL's own app data: logs, single-instance lock, WebView2 profile. Never the user's results
; (those are next to their input files or wherever they saved downloads).
Type: filesandordirs; Name: "{localappdata}\COAL"
; Anything the app left inside its own folder.
Type: filesandordirs; Name: "{app}"
; CHANGED(WIN-W6): new. The local app's run folders (results shown in the window; downloads you
; saved elsewhere are untouched). See overnight/windows-app/CHANGES.md
Type: filesandordirs; Name: "{%TEMP}\coal_app_runs"

[Code]
// Microsoft's documented way to detect the Evergreen WebView2 Runtime: the "pv" value under
// EdgeUpdate\Clients\{F3017226-FE2A-4295-8BDF-00C3A9A7E4C5}, machine-wide (WOW6432Node) or
// per user. COAL still installs without it, but warns, because the window needs it.
const
  WebView2Key = 'Microsoft\EdgeUpdate\Clients\{F3017226-FE2A-4295-8BDF-00C3A9A7E4C5}';

function WebView2Installed(): Boolean;
var
  Version: String;
begin
  Result :=
    (RegQueryStringValue(HKLM, 'SOFTWARE\WOW6432Node\' + WebView2Key, 'pv', Version) and (Version <> '') and (Version <> '0.0.0.0')) or
    (RegQueryStringValue(HKLM, 'SOFTWARE\' + WebView2Key, 'pv', Version) and (Version <> '') and (Version <> '0.0.0.0')) or
    (RegQueryStringValue(HKCU, 'Software\' + WebView2Key, 'pv', Version) and (Version <> '') and (Version <> '0.0.0.0'));
end;

function InitializeSetup(): Boolean;
begin
  Result := True;
  if not WebView2Installed() then
    MsgBox('COAL shows its window with the Microsoft Edge WebView2 Runtime, which this computer '
      + 'does not seem to have (it is part of Windows 11 and of Edge updates on Windows 10).'
      + #13#10#13#10 + 'COAL will still install. If its window does not open, COAL offers to open in '
      + 'your web browser instead (still only on this computer), or you can install the '
      + '"Evergreen Bootstrapper" from Microsoft''s WebView2 page.', mbInformation, MB_OK);
end;

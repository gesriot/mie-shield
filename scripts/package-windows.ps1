param(
    [string]$PythonVersion = "3.14",
    [string]$ProductName = "MieShield",
    [string]$BuildEnv = ".build-venv/windows",
    [string]$OutputDir = "dist/windows"
)

$ErrorActionPreference = "Stop"
$RootDir = Split-Path -Parent $PSScriptRoot
Set-Location $RootDir

$env:UV_CACHE_DIR = if ($env:UV_CACHE_DIR) { $env:UV_CACHE_DIR } else { Join-Path $RootDir ".uv-cache" }
$env:NUITKA_CACHE_DIR = if ($env:NUITKA_CACHE_DIR) { $env:NUITKA_CACHE_DIR } else { Join-Path $RootDir ".nuitka-cache" }
$env:MPLCONFIGDIR = if ($env:MPLCONFIGDIR) { $env:MPLCONFIGDIR } else { Join-Path $RootDir ".mplconfig" }

$BuildEnvPath = Join-Path $RootDir $BuildEnv
$OutputPath = Join-Path $RootDir $OutputDir

Write-Host "==> Creating build environment: $BuildEnvPath"
uv venv --python $PythonVersion $BuildEnvPath

$Version = & "$BuildEnvPath/Scripts/python.exe" -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])"
$IconPath = Join-Path $RootDir "dist/icons/$ProductName.ico"
$ZipPath = Join-Path $RootDir "dist/$ProductName-$Version-windows-x64.zip"

Write-Host "==> Installing runtime dependencies into build environment"
$env:UV_PROJECT_ENVIRONMENT = $BuildEnvPath
uv sync --locked

Write-Host "==> Installing build-only dependencies"
uv pip install --python "$BuildEnvPath/Scripts/python.exe" "nuitka>=4,<5" ordered-set zstandard
Remove-Item Env:\UV_PROJECT_ENVIRONMENT

if (Test-Path "icon.png") {
    Write-Host "==> Generating Windows icon from icon.png"
    & "$BuildEnvPath/Scripts/python.exe" scripts/make_icon.py --input icon.png --ico $IconPath
}

Write-Host "==> Removing previous Windows build output"
Remove-Item -Recurse -Force "$OutputPath/mie_shield.dist", "$OutputPath/mie_shield.build", $ZipPath -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force $OutputPath | Out-Null

$NuitkaArgs = @(
    "--mode=standalone",
    "--windows-console-mode=disable",
    "--enable-plugin=pyside6",
    "--include-package=PyMieScatt",
    "--lto=yes",
    "--python-flag=-O",
    "--python-flag=no_docstrings",
    "--output-filename=$ProductName.exe",
    "--output-dir=$OutputPath"
)

if (Test-Path $IconPath) {
    $NuitkaArgs += "--windows-icon-from-ico=$IconPath"
}

Write-Host "==> Building Windows application with Nuitka"
& "$BuildEnvPath/Scripts/python.exe" -m nuitka @NuitkaArgs mie_shield.py

Write-Host "==> Creating ZIP artifact"
Compress-Archive -Path "$OutputPath/mie_shield.dist/*" -DestinationPath $ZipPath -Force

Write-Host "==> Done"
Write-Host "Dist: $OutputPath/mie_shield.dist"
Write-Host "ZIP:  $ZipPath"

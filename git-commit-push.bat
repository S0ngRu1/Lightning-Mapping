@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion

:: 检查是否在 Git 仓库中
git rev-parse --git-dir >nul 2>nul
if errorlevel 1 (
    echo ❌ 错误：当前目录不是 Git 仓库！
    echo 请在项目根目录（包含 .git 文件夹）运行此脚本。
    echo 例如：cd /d D:\your-project && git-commit-push.bat "你的提交信息"
    exit /b 1
)

:: 获取提交信息
set "commit_msg=Update"
if not "%～1"=="" (
    set "commit_msg=%*"
)

echo.
echo 📥 正在添加所有更改...
git add . >nul 2>nul

:: 检查是否有更改
git diff --cached --quiet
if errorlevel 1 (
    echo 📝 正在提交: %commit_msg%
    git commit -m "%commit_msg%" --quiet
    if errorlevel 1 (
        echo ❌ 提交失败！请检查：
        echo   1. 未添加任何文件（git add .）
        echo   2. 提交信息是否为空
        exit /b 1
    )
) else (
    echo ⚠️ 无文件更改，跳过提交。
    exit /b 0
)

echo.
echo 🚀 正在推送...
git push --quiet

:: 修复关键：正确检查 errorlevel
if errorlevel 1 (
    echo ❌ 推送失败！错误详情：
    echo   1. 网络问题（SSL 错误：SSL_ERROR_SYSCALL）
    echo   2. 请检查：
    echo      - 代理设置（如 Clash/v2rayN 是否开启局域网共享）
    echo      - Git 证书（重新安装 Git for Windows 时勾选 "Use OpenSSL"）
    echo      - 远程地址：git remote -v
    echo      - 尝试切换为 SSH 地址：git remote set-url origin git@github.com:username/repo.git
    exit /b 1
) else (
    echo ✅ 提交并推送成功！
    exit /b 0
)
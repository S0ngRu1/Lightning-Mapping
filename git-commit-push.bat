@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion

:: 检查是否在 Git 仓库中
git rev-parse --git-dir >nul 2>nul
if errorlevel 1 (
    echo ❌ 错误：当前目录不是 Git 仓库！
    echo 请在项目根目录（包含 .git 文件夹）运行此脚本。
    exit /b 1
)

:: 获取提交信息（支持带空格的参数）
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
) else (
    echo ⚠️ 无文件更改，跳过提交。
    exit /b 0
)

echo.
echo 🚀 正在推送...
git push --quiet

if errorlevel 0 (
    echo ✅ 提交并推送成功！
    exit /b 0
) else (
    echo ❌ 推送失败！请检查：
    echo   1. 网络连接
    echo   2. 远程仓库地址 (git remote -v)
    echo   3. 本地分支与远程分支是否匹配
    exit /b 1
)
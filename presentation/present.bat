@echo off
cd /d "%~dp0"
call "..\.venv\Scripts\activate.bat"
manim-slides present --show-info-window Preliminaries Heart Isogenies

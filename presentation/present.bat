@echo off
cd /d "%~dp0"
call "..\.venv\Scripts\activate.bat"
manim-slides present --start-paused Preliminaries Isogenies

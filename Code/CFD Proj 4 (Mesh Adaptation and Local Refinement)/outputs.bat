@echo off
set /p mesh=Mesh: 

@REM rem compile solver
g++ -o main.exe 1st_order.cpp
g++ -o finres.exe finres.cpp

python python/plotgri.py %mesh%
python python/plotmach.py %mesh% mach
python python/calc.py %mesh%
pause
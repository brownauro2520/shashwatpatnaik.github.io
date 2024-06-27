@echo off
set /p num_runs=Enter the number of times to run the commands:

@REM rem compile solver
g++ -o main.exe 1st_order.cpp
g++ -o finres.exe finres.cpp

for /l %%i in (3,1,%num_runs%) do (
    @REM rem make matrices for coarse mesh
    python python/makemats.py %%i coarse

    @REM rem converge solver
    main.exe %%i inf coarse

    @REM rem inject the coarse state into a finer mesh, save matrices and create finer grid
    python python/inject.py %%i

    @REM rem compute the fine space injected residual RUHh
    finres.exe %%i

    @REM rem converge the fine space to steady state
    main.exe %%i load fine

    @REM compute output jacobian
    python python/ojac.py

    @REM compute adjoint using Nathans code
    @REM python python/natadj.py

    @REM flag % of elements and do refinement
    python python/localref.py %%i
)

pause
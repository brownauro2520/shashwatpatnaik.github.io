@echo off
set CPP_FILE=rk2.exe
set CPP_ARG1=0
set CPP_ARG2=0.2
set CPP_ARG3=1
set CPP_ARG4=0.1
set CPP_ARG5=2
set CPP_ARG6=0.1
set CPP_ARG7=3
set CPP_ARG8=0.1



set PYTHON_FILE=test.py
set PYTHON_ARG1=0
set PYTHON_ARG2=1
set PYTHON_ARG3=2


rem Run the C++ program
"%CPP_FILE%" "%CPP_ARG1%" "%CPP_ARG2%" 

rem Run the Python program
python "%PYTHON_FILE%" "%PYTHON_ARG1%" 

rem Run the C++ program
"%CPP_FILE%" "%CPP_ARG3%" "%CPP_ARG4%" 

rem Run the Python program
python "%PYTHON_FILE%" "%PYTHON_ARG2%"

rem Run the C++ program
"%CPP_FILE%" "%CPP_ARG5%" "%CPP_ARG6%" 

rem Run the Python program
python "%PYTHON_FILE%" "%PYTHON_ARG3%"

rem Run the C++ program
"%CPP_FILE%" "%CPP_ARG7%" "%CPP_ARG8%" 

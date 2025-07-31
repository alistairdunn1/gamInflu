@echo off
rem =====================================
rem Build and Check gamInflu R Package (Windows)
rem =====================================

rem Check if R is in PATH
where R >nul 2>nul
if errorlevel 1 (
  echo R is not in your PATH. Please install R and add it to your PATH.
  exit /b 1
)

rem Generate documentation with roxygen2
echo Generating documentation...
call R --vanilla < run-roxygen.R
if errorlevel 1 exit /b 1

rem Build the package
echo Building package...
call R CMD build --force gamInflu
if errorlevel 1 exit /b 1

rem Install the package
echo Installing package...
call R CMD INSTALL gamInflu
if errorlevel 1 exit /b 1

rem Check the package
echo Checking package...
call R CMD check gamInflu
if errorlevel 1 exit /b 1

echo.
echo =====================================
echo Build completed successfully!
echo =====================================

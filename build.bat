rem Build R library and version
R --vanilla < run-roxygen.R

rem Build library
rem R CMD build --force GAMInflu
rem R CMD INSTALL --build GAMInflu
rem R CMD check GAMInflu

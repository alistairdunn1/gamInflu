rem Build R library and version
R --vanilla < run-roxygen.R

rem Build library
R CMD build --force GAMInflu
R CMD INSTALL --build GAMInflu
R CMD check GAMInflu

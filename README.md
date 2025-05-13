# GAMInflu
GAM Influence plots

# R code for removing and reinastalling after a new build
detach("package:GAMInflu", unload = TRUE)
remove.packages("GAMInflu")
install.packages("C:/Users/alist/OneDrive/Projects/Software/GAMInflu/GAMInflu_0.1.0.tar.gz", clean = TRUE)
library(GAMInflu)
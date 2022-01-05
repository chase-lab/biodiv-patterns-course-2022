## Morning session

Metric introduction (?)
GitHub introduction

## Afternoon session  

Please, BEFORE THE CLASS, run the following R code to install packages in preparation of our exercise with MoBsim:
```
install.packages(c('shiny', 'shinyBS','pals','shinyjs', 'devtools'))
devtools::install_github('albansagouis/mobsim@master')
```
Please, DO NOT install `mobsim` through `install.packages("mobsim")`.

And make sure this works:
```
shiny::runGitHub("albansagouis/mobsim_app", ref = "master")
```
If the shiny app does not appear or if it crashes, please send me an email: alban.sagouis@idiv.de


## Morning session

Class beginning time: 09:30  

### Introduction (to an introduction) to GitHub
Please, BEFORE THE CLASS:
 + Create an account on [github.com][https://www.github.com]
 + Install GitHub Desktop.  
   - For Windows and MacOS: https://desktop.github.com/  
   - For Linux: https://www.vitainbeta.org/2020/06/24/installing-github-desktop-linux/


### MoBsim
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


## Afternoon session

Description of the projects to prepare for Friday afternoon.  
Group assignments, and help getting started (we're there to help you)

Some code about data extraction in R can be found in this page: https://github.com/chase-lab/biodiv-patterns-course-2021/tree/main/week%202/3%20-%20Wednesday

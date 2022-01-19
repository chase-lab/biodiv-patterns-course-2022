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

What can we learn from studying different biodiversity metrics? 
-	Richness (S or alpha): number of species, but will be heavily influences by rare species with an abundance of 1.  
-	Abundance (N): are there differences in number of individuals?
-	Rarefied richness: can the difference in richness simply be explained by differences in abundance? 
-	PIE – is there a difference in evenness between treatments? (note that Shannon is always somewhere in between alpha and PIE, and may therefore not be super interesting. 
-	Gamma diversity in comparison to alpha diversity: is the difference in richness scale dependent (is the interpretation for Gamma different than for Alpha?). Could this indicate homogenization of the community? 
-	Beta diversity: is treatment X more homogeneous than treatment Y? is there a difference in species composition between treatments? 



Group assignment 

Find 5 papers that compared different treatments of some ecological factor (land use, different time points, elevation, pollution, etc) in terms of species richness or other metrics and calculate our other favorite metrics: Rarefied richness, PIE, beta (if you really want to diver into that rabbit hole), etc. How often do you find different interpretations than the original authors?  
You will need raw data that provide for each species the abundance in each plot (not aggregated over plots or treatments). Not all papers provide this. Look for papers published in the last 3-4 years. 

1)	Find some papers that study species richness differences 
2)	Download the underlying raw data 
3)	Wrangle the data to put it in the right format 
4)	Check the result of the original papers
5)	Calculate the metrics you're interested in (alpha, beta, gamma, rarefied richness, PIE, Shannon etc)
6)	Compare to the results of the authors. How often do you find differences? Which important patterns/ differences have the original authors missed (if any)?
 

2 paths: 
Find publication (Google scholar, Web of Science (university access)) -> data repo 
Repositories (Dryad, NERC, EDI, LTER) -> find paper that belongs to it. 

Some code about data extraction in R can be found in this page:  
https://github.com/chase-lab/biodiv-patterns-course-2021/tree/main/week%202/3%20-%20Wednesday



What we expect for the presentation: 

1) question
2) Methods (explanation of datasets, how did you analyse them? which metrics did you extract? etc.) 
3) Results with graphs 
4) Conclusions.







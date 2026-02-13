# BVSMed_GLM

You can check the results of Scenario 4 of simulation in this paper through the code posted here.

When you download the folder, a folder called BVSMed_GLM-main is created and Rstudio is executed under this folder.
This means that the working directory of the Rstudio should be set as follows: "~/BVSMed_GLM-main"


The data imported from the code was created through the data folder's data.R code.
The eta value used within the code was determined based on the rules presented in the paper based on the results of the find_eta.R code.
The description of the arguments of the function is written in the code simulation.R.

In the paper, we focus on a binary outcome with a logit link. 
Accordingly, in the main code, link = 1 specifies the logit link for binary outcomes. 
To use count data with a log link (Poisson regression), set link = 2.

# This script has one purpose; to render the Rmd file into the correct directory and in the correct format.

library(rmarkdown)
library(here)
render(here("src", "ExoMeth_Manuscript.Rmd"), 
       output_dir = here("docs"))

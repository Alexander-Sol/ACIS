# Quick Convert, just change the doc name
library(knitr)

doc.name <- "AutoClustR_Outline"

# doc.name <- readline(prompt = "Enter the path to the .docx file to be converted (do not include the file extension)") 

if(!is.na(file.info(paste0(doc.name, ".markdown")))) {
  file.remove(paste0(doc.name, ".markdown"))
}

knitr::pandoc(input = paste0(doc.name, ".docx"),
              format = "markdown")



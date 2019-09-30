library(holepunch)

options(usethis.full_name = "Shea Connell", 
        usethis.description = list(`Authors@R` = 'person("Shea", "Connell", email = "sheaconnell@gmail.com", role = c("aut", "cre"), 
                                   comment = c(ORCID = "0000-0002-2850-3908"))', 
                                   License = "MIT + file LICENSE", 
                                   Version = "0.0.1"))

write_compendium_description(package = "The Framework: ExoMeth", 
                             description = "This contains all of the data and code neccessary to reproduce the ExoMeth manuscript analyses, 
                             figures and the majority of the formatted manuscript")

write_dockerfile(maintainer = "Shea Connell")

generate_badge()

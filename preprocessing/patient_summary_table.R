# patient summary table

# median age
median(mutation_covariate_data$age)

# median age per cancer type
median(mutation_covariate_data$age[mutation_covariate_data$Dx=="AML"])
median(mutation_covariate_data$age[mutation_covariate_data$Dx=="MDS"])
median(mutation_covariate_data$age[mutation_covariate_data$Dx=="CMML"])
median(mutation_covariate_data$age[mutation_covariate_data$Dx=="MPN"])

# range age
range(mutation_covariate_data$age)

# median age per cancer type
range(mutation_covariate_data$age[mutation_covariate_data$Dx=="AML"])
range(mutation_covariate_data$age[mutation_covariate_data$Dx=="MDS"])
range(mutation_covariate_data$age[mutation_covariate_data$Dx=="CMML"])
range(mutation_covariate_data$age[mutation_covariate_data$Dx=="MPN"])


# median age
sum(mutation_covariate_data$Gender=="Male")

# median age per cancer type
sum(c(mutation_covariate_data$Gender=="Male")[mutation_covariate_data$Dx=="AML"])
sum(c(mutation_covariate_data$Gender=="Male")[mutation_covariate_data$Dx=="MDS"])
sum(c(mutation_covariate_data$Gender=="Male")[mutation_covariate_data$Dx=="CMML"])
sum(c(mutation_covariate_data$Gender=="Male")[mutation_covariate_data$Dx=="MPN"])


# median age
sum(mutation_covariate_data$Gender=="Female")

# median age per cancer type
sum(c(mutation_covariate_data$Gender=="Female")[mutation_covariate_data$Dx=="AML"])
sum(c(mutation_covariate_data$Gender=="Female")[mutation_covariate_data$Dx=="MDS"])
sum(c(mutation_covariate_data$Gender=="Female")[mutation_covariate_data$Dx=="CMML"])
sum(c(mutation_covariate_data$Gender=="Female")[mutation_covariate_data$Dx=="MPN"])



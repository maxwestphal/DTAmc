stats <- data2stats(data, covariance = T, regu=T)



set.seed(1337)
data <- sample_data(m=10)

dta(data, test="univariate")
dta(data, test="bonf")
dta(data, regu=F, test="maxt")
dta(data, regu=T, test="maxt")
dta(data, test="maxt", useSEPM=T)

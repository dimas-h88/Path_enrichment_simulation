library(Matrix)
sigma <- matrix(runif(10000^2, 0, 0.1), nrow = 10000)
  
for (i in 1:10000) {
  sigma[i, i] <- 0.25
}
  
sigma <- nearPD(sigma)$mat
sigma <- as.matrix(sigma)  # Convert to a regular matrix
  
write.csv(sigma, file = 'sigma_matrix_0.1.csv')

sigma <- matrix(runif(10000^2, 0, 0.01), nrow = 10000)

for (i in 1:10000) {
  sigma[i, i] <- 0.25
}

sigma <- nearPD(sigma)$mat
sigma <- as.matrix(sigma)  # Convert to a regular matrix

write.csv(sigma, file = 'sigma_matrix_0.01.csv')


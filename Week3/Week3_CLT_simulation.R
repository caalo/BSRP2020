N <- 100000 # population size
# create our population 
pop_dist1 <- rnorm(N, mean = 50, sd = 10) # from a normal dist
pop_dist2 <- rnorm(N, mean = 20, sd = 10) # from a normal dist
#pop_dist <- runif(N, min = 0, max = 100) # from a uniform dist
pop_dist <- rpois(N, lambda = c(25, 50, 75))
hist(pop_dist)
mean(pop_dist)
sd(pop_dist)

n <- 500 # sample size
n_samples <- 10000 # how many times are we sampling?


sample_means <- numeric(length = n_samples) # here, we're going to store the means of each sample

# sample the number of times as n_samples
for(i_sample in 1:n_samples){
  sample_i = sample(pop_dist, size = n) # sample size = n
  sample_means[i_sample] <- mean(sample_i) # store the mean of the sample into sample_means
}

hist(sample_means)

mean(sample_means)
sd(sample_means)


#STATISTICS AND HYPOTHESIS TESTING DEMO

#Motivation: What is a population? What is a sample? Why might a sample not be the "same" as the population?

#suppose that you are doing a study on sepal length measurements of iris plants in Boston
#what is a population distribution?

iris_population <- iris
dim(iris_population)
hist(iris_population$Sepal.Length, 100)

#as scientists, what are some ways to describe/summarize a population's distribution?

mean(iris_population$Sepal.Length)
sd(iris_population$Sepal.Length)
median(iris_population$Sepal.Length)

#what is a sample distribution?

#how do we get a sample distribution from a population distribution?

iris_sample = iris_population[sample(1:nrow(iris_population), 20) ,]

#what are some ways to summarize/describe a sample distribution? is it a random variable?
##what are we trying to understand with this value?

mean(iris_sample$Sepal.Length)
iris_sample = iris_population[sample(1:nrow(iris_population), 20) ,]
mean(iris_sample$Sepal.Length)
iris_sample = iris_population[sample(1:nrow(iris_population), 20) ,]
mean(iris_sample$Sepal.Length)
iris_sample = iris_population[sample(1:nrow(iris_population), 20) ,]
mean(iris_sample$Sepal.Length)

#what are some factors that can change the sample mean?
iris_sample = iris_population[sample(1:nrow(iris_population), 50) ,]
mean(iris_sample$Sepal.Length)


#Let's do this a bunch of times: want to see whether the sample mean is consistently similar to the population mean
#
set.seed(0123)
n_expts = 10000
sample_size = 20
sample_means = numeric()
for(i in 1:n_expts) {
  iris_sample = iris_population[sample(1:nrow(iris_population), sample_size) ,]
  sample_means <- c(sample_means, mean(iris_sample$Sepal.Length))
}

iris_sample = iris_population[sample(1:nrow(iris_population), sample_size) ,]

hist(sample_means, 100) 
#Does the peak of the sample mean distribution correspond to the population mean?
abline(v = mean(iris_sample$Sepal.Length), col = "red", lwd = 3)
abline(v = mean(iris_population$Sepal.Length), col = "blue", lwd = 3)
abline(v = mean(sample_means), col = "green", lwd = 3)

mean(iris_population$Sepal.Length)
mean(sample_means)


#in reality, we only get the sample ONCE (usually), 
#and we don't know what the population mean is. 
#so, we can't do the following: 
#compare the sample mean to population mean 
#compare the sample mean to the distribution of the sample means. 

#it seems that we can't easily figure out what the population mean is going to be, 
#given the sample mean.

#A fact from statistics: given a population mean, with a sample size,
#we can predict what the sample mean distribution looks like. 
#For instance, a population mean of 0, take a sample size of 20:
t_distribution <- rt(n = 10000, df=20)
hist(t_distribution, 100)
#then, suppose we observe a mean of 2 in our sample:
abline(v = 2, col = "red", lwd = 3)
#what is the probability of observing 2 or any value more extreme?
mean(t_distribution > 2)
#not very likely (< .05), so the population mean is not very likely to be 0.

#We have, instead, worked backwards: instead of staring at the sample mean and trying
#to deduce the population mean, we assumed what the population mean is a particular value,
#deduced its properties, and see whether the sample mean fit the property. 

#This is the concept of hypothesis testing. 
#Have a hypothesis what the population (mean) is. 
#The sample mean follows a distribution: the null distribution. 
#Take a sample and look at it: where does it fall in the in the null distribution?
#If it is not a very likely value in the null distribution, then it is likely that
#the population mean is not that hypothesis.

#Back to the iris example.
#Suppose null hypothesis: population mean is 5. 
null_hypothesis <- 5
t_distribution <- rt(10000, df = nrow(iris_sample) - 1)
hist(t_distribution, 100)
t_statistic <- (mean(iris_sample$Sepal.Length) - null_hypothesis) / sqrt(var(iris_sample$Sepal.Length) / nrow(iris_sample))
abline(v = t_statistic, col = "red", lwd = 3)
mean(t_distribution > t_statistic)
t.test(iris_sample$Sepal.Length, mu = 5)
#try another null hypothesis!

#Now, what if the null hypothesis is the population mean of 5.84?
null_hypothesis <- 5.84
t_statistic <- (mean(iris_sample$Sepal.Length) - null_hypothesis) / sqrt(var(iris_sample$Sepal.Length) / nrow(iris_sample))
t_distribution <- rt(10000,df = nrow(iris_sample) - 1)
hist(t_distribution, 100)
abline(v = t_statistic, col = "red", lwd = 3)
mean(t_distribution > t_statistic)

#Is it possible that by chance we observe something that will happen less than 5% of the time?

#############

#IN-CLASS EXERCISE:

#we came up with the following scientific hypothesis: there is a difference between sepal length between
#setosa and virginica. 
library(tidyverse)
iris_setosa_population <- iris_population %>% filter(Species == "setosa")
iris_virginica_population <- iris_population %>% filter(Species == "virginica")

population_difference <- mean(iris_setosa_population$Sepal.Length) -
  mean(iris_virginica_population$Sepal.Length)

#Take a sample of 20 from iris_setosa_population and iris_virginica_population. 
#What is the difference in means?

#The difference in population means actually follows a T-distribution also.

#Based on our original scientific hypothesis, what should be the null hypothesis be?

#Here is the way to generate the null hypothesis and t_statistic 
#(with variables iris_setosa_sample, iris_virginica_sample):
t_distribution <- rt(n = 10000, df = 76)
t_statistic <- (mean(iris_setosa_sample$Sepal.Length) -  mean(iris_virginica_sample$Sepal.Length)) / 
  sqrt(var(iris_setosa_sample$Sepal.Length) / 20 + var(iris_virginica_sample$Sepal.Length) / 20)

#what is the p-value?
#verify with t.test()!

#Look back at Q4 of the homework. Can you do it now?
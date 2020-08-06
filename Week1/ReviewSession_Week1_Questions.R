# QUESTIONS #

#### QUESTION 1: ####
# Which of the following lines of code gives me the 3rd column of the iris dataset?

# A:
iris[3]
# B:
iris(3)
# C:
iris[3,]
# D:
iris[,3]
# E:
# NONE OF THE ABOVE

#### QUESTION 2: ####
# Which of the following loops correctly prints the numbers 1-10?

# A:
for(a in c(1,10)){
  print(a)
}
# B:
b=1
while(b <= 10){
  print(b)
  b=b+1
}
# C:

add2 <- function(x){
  x + 2
}
sapply(1:10, add2)

# D:
print(1:10)

## Q



#### QUESTION 3: ####
# How many times does the following code chunk print "*"?
i = 1
while(i < 10){
  print("*")
  if(i == 5){
    i = i + 2
  }else{
    i = i + 1
  }
}

# A: 1
# B: 5
# C: 9
# D: 10
# E: None of the above
#### QUESTION 4: ####

# Here, I generate 100 random numbers from a normal distribution (mean = 0, sd = 1)
random_numbers = rnorm(100)
# How can I make a new vector containing the numbers less than 0?

# A:
new_vector <- c(random_numbers < 0)

which(random_numbers < 0)

# C:
new_vector <- random_numbers[random_numbers < 0]



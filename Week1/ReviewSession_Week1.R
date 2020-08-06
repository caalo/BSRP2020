
### SUBSET A VECTOR ####
myvector <- c(1,2,3,4,5,1,2,3,4,5)
#             1 2 3 4 5 6 7 8 9 10

# get values less than 3
which(myvector < 3)
myvector[which(myvector < 3)]
myvector[myvector < 3]

myvector[c(1,2,6,7)]

### DOING IT THE LONG WAY ###
newvector <- numeric()
index_newvector <- 1

for(index_myvector in 1:length(myvector)){
  if(myvector[i] < 3){
    newvector[j] <- myvector[i]
    j = j+1
  }
}

### DATA TYPES ###
# Vectors
# Initializing a vector vs converting data into a vector
my_numeric_vector <- numeric(10)
my_numeric_vector_2 <- as.numeric(10)
my_numeric_vector_2.5 <- as.numeric("10")
my_numeric_vector_3 <- as.numeric(1:10)
is.numeric(my_numeric_vector_3)

my_character_vector <- character(10)
my_character_vector <- character("stuff")
my_character_vector_2 <- as.character("stuff")
as.character(10)

my_vector <- c(1,2,3,4,5,6,7,8,9,10)
my_named_vector <- c("a" = 1, "b" = 2, "c" = 3)
class(my_vector)
class(my_named_vector)

my_vector <- as.character(my_vector)
class(my_vector)

my_vector <- as.character(my_vector)
class(my_vector)

my_vector2 <- c("A", "T", "G", "C")
class(my_vector2)

as.character(my_vector2)
as.numeric(my_vector2)
class(my_vector2)

my_list <- list(
  a = c(1,2,3),
  b = c(4,5,6),
  c = c(7,8,9),
  d = c(10,11,12)
)

deo_list <- list(Gisselle = c(1,2,3,4), 
                 Sid = c(0,3), 
                 Chris = c(2,4,6,8))

my_list[[1]]
my_list[["a"]]
antonio = my_list[1]
my_list[["d"]][2]

new_list <- list(
  list1 = my_list,
  list2 = deo_list
)

new_list[["list1"]][["d"]][2]
new_list$list1$d[2]

lapply(new_list, )


### LIST TO DATAFRAME? ###
as.data.frame(my_list)

### DATAFRAME TO LIST? ###
as.list(iris)

### STRING SPLITTING ####
sentence <- "Hello world, it's me, Bob"
two_sentences <- c("Hello world, it's me, Bob.", "Did you miss me?")

# strsplit applies itself to each element of a list/vector
split_sentence <- strsplit(sentence, split = " ")
split_2sentences <- strsplit(two_sentences, split = " ")

# strsplit can also split a single word into individual characters, if split is just ""
split_word <- strsplit("sentence", split = "")

# We can access the elements of the list that strsplit produces
split_sentence[1] # first element of the list
split_sentence[[1]] # also the first element of the list
split_sentence[[1]][1] # first element of the first element of the list

### SUBSETTING A DATAFRAME

# Here, we're using the in-built iris dataset. You don't need to load this in, it already exists in R.
# There are many ways of subsetting a dataframe
# The first is the matrix-like way: we can access rows and columns, using [,]
# whatever comes before the , is the row numbers, and whatever comes after is the columns.
iris[1,] # first ROW of the iris data set. Includes ALL columns
iris[,1] # first COLUMN of the iris data set. Includes ALL rows

iris[c(1,2),c(3,4)] # contents of the 1st and 2nd rows, and 3rd and 4th columns
iris[1:5,1:5] # contents of first 5 rows, and first 5 columns

iris[-(1:5),] # everything EXCEPT the first 5 rows

# We can also access specific column names - as long as they are quoted!
iris[,"Sepal.Width"]

# The second way of subsetting is the list-like way: we can access specific columns using [[]]
iris[["Sepal.Width"]]
# Note: when we use the list-like way, we can't specify rows. This only lets us use columns.

# The third way of subsetting is the $. This lets us access a colunn WITHOUT quoting it.
iris$Sepal.Width

# We can combine the different subsetting methods to give us various slices of data
iris[iris$Sepal.Width < 2.5, "Species"] # Take rows where Sepal.Width < 2.5, and take the "Species" column
iris[iris$Sepal.Width < 2.5, c("Species", "Sepal.Width")] # Same thing, but take both the "Species" and "Sepal,Width" columns

iris[["Species"]][iris$Sepal.Width < 2.5] # same as line 122 - but we are subsetting from the extracted column, which becomes a vector
iris$Species[iris$Sepal.Width < 2.5]

### FUNCTIONS ###

multiply_by_3 <- function(input){
  result = input * 3
  return(result)
}

get_nth_word_of_sentence <- function(string, n = 3){
  words_list = strsplit(string, split = " ")
  
  danger_length = length(words_list[[1]])+1
  if(n >= danger_length){
    stop("Error - you went too far. Turn back now!")
  }else{
    word = words_list[[1]][n]
    return(word)
  }
}

sentence <- "Hello world, it's me, Bob"
get_nth_word_of_sentence(sentence, n = 3)

list_of_sentences <- list(
  a = "Hello world",
  b = "It's a-me",
  c = "Mario!!!!",
  d = "The crazy plumber"
)

lapply(list_of_sentences, get_nth_word_of_sentence)



### LOOPS ###

sofias_vector <- c(0)

add_to_sofias_vector <- function(new_element, sofias_vector){
  sofias_vector <- c(sofias_vector, new_element)
  return(sofias_vector)
}

sofias_vector <- c(0)

sofias_vector <- lapply(1:10, add_to_sofias_vector, sofias_vector)

i = 1
while(i < 10){
  print(i)
  i = i+1
}

sapply(1:10, function(x){
  return(x * 2)
})

lapply(1:10, print)

sapply(1:10, multiply_by_3)

#### FUNCTION + LOOP ####
triangle_print <- function(height = 3){
  for(i in 1:height){
    stars_i = rep("*", i)
    row_i = paste0(stars_i, collapse = "")
    print(row_i) # cat is like print
  }
  for(i in (height-1):1){
    stars_i = rep("*", i)
    row_i = paste0(stars_i, collapse = "")
    print(row_i) # cat is like print
  }
}

### DNA TRANSCRIPTION ###

transcribe <- function(dna_string){
  # split the string up into individual letters
  
  # Write a conditional statement - if each letter is an A, T, C, or G,
  # return a new letter that corresponds to the transcribed RNA nucleotide
  # e.g. U, A, G, of C.
  
  # PASTE all of the letters together to form a new string called rna_string
  
  return(rna_string)
}















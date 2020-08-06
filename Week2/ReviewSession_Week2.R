### WEEK 2 Review ###
library(tidyverse)

#### select ####
iris %>%
  select(Species)

iris %>%
  select(starts_with("Petal"))

iris %>%
  select(-starts_with("Sepal"))

#### mutate ####
iris %>% 
  mutate(has_setosa = Species == "setosa")

# Same as:
cbind(iris, has_setosa = iris$Species == "setosa")
# or:
iris$has_setosa <- iris$Species == "setosa"

#### unite ####
CCLE_metadata <- read.csv("sample_info.csv") # read in the CCLE metadata

CCLE_metadata %>%
  unite("new_names", stripped_cell_line_name, lineage) %>%
  glimpse()

#### separate ####
CCLE_metadata %>%
  separate(CCLE_Name, sep = "_", into = c("CellLine", "TYPE"), extra = "merge") %>%
  glimpse()

#### rename ####
CCLE_metadata %>%
  rename(CCLE_ID = CCLE_Name) %>%
  glimpse()

#### recode ####
iris %>%
  recode(Species, setosa = "setty")

#### gather ####
iris_gathered <- iris %>%
  rowid_to_column(var = "id") %>% # add row ids to gather (each row must have a unique identifier for gather to work)
  gather(measure_type, measure_value, -id, -Species)

#### spread ####
iris_spread <- iris_gathered %>%
  spread(measure_type, measure_value) %>%
  select(-id) # drop the id column we added from the rowids

#### group_by ####
iris %>%
  group_by(Species)

#### summarize ####
iris %>%
  group_by(Species) %>%
  summarize(mean(Petal.Length))

#### merge - actually part of base R####
iris_descriptions <- data.frame(
  Species = c("setosa", "versicolor", "virginica"),
  description = c("tall, dark", "short, bright", "has nothing to do with virginia")
)

# NOTE: BOTH DATAFRAMES MUST HAVE 1 OR MORE COLUMN NAMES IN COMMON
iris_merged <- merge(iris_descriptions, iris)

#### match / %in% - actually part of base R ####
match(1.6, iris$Petal.Length) # gives the index of the first occurrence of a value
1.6 %in% iris$Petal.Length # checks to see if the value is in a vector
iris$Petal.Length %in% 1.6 # checks to see if the vector contains the value.

which(iris$Petal.Length %in% 1.6) # which values in the vector are equal to 1.6
which(iris$Petal.Length %in% c(1.4, 1.6)) # which values in the vector are equal to 1.6 or 1.4

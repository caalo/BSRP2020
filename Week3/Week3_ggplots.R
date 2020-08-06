# Review Session - Week 3
library(tidyverse)

#### EXAMPLE 1 ###
# Wrangling the iris data to plot petal length setosa vs versicolor
iris_toPlot <- iris %>%
  rowid_to_column() %>%
  select(Petal.Length, Species, rowid) %>%
  spread(Species, Petal.Length)

#### EXAMPLE 2 ####
# Adding shapes/colors to ggplot
ggplot(data = iris, aes(x = Petal.Length, y = Petal.Width, color = Species)) +
  geom_point()

ggplot(data = iris, aes(x = Petal.Length, y = Petal.Width, shape = Species)) +
  geom_point()

# Adding a single color/shape
ggplot(data = iris, aes(x = Petal.Length, y = Petal.Width, shape = Species)) +
  geom_point(shape = 3)

ggplot(data = iris, aes(x = Petal.Length, y = Petal.Width, shape = Species)) +
  geom_point(color = "purple", size = 2)

#### EXAMPLE 3 ####
iris_plot <- ggplot(data = iris, aes(x = Petal.Length, y = Petal.Width, shape = Species, color = Species)) +
  geom_point()


ggplot(data = iris, aes())
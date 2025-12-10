Pumpkins <- read.csv("pumpkins_16.csv")
library(tidyverse)
the_heaviest <- Pumpkins %>%
  filter(weight_lbs == max(weight_lbs, na.rm = TRUE)) %>%
  select(variety, city, country, id, weight_lbs)
the_heaviest

weight <- function(pounds_to_kilogram){
  pounds_to_kilogram*0.45
}
Pumpkins$weight_kg <- weight(Pumpkins$weight_lbs)
head(Pumpkins)

#2 Identify the heaviest pumpkin grown in any of the competitions – what variety was it? Where was it from and when was it grown?

library(tidyverse)
the_heaviest <- Pumpkins %>%
  filter(weight_lbs == max(weight_lbs, na.rm = TRUE)) %>%
  select(variety, city, country, id, weight_lbs)
the_heaviest

#3 Write a function to change the weight in pounds (lbs) to kilograms (kg) and use this function to create a new column in your data set called weight_kg.
weight <- function(pounds_to_kilogram){
  pounds_to_kilogram*0.45
}
Pumpkins$weight_kg <- weight(Pumpkins$weight_lbs)
head(Pumpkins)

#4 Create another new column in your data set called weight_class. 
Pumpkins$weight_class = ifelse(Pumpkins$weight_lbs < 1000, "light", ifelse(Pumpkins$weight_lbs < 1500, "medium","heavy"))
head(Pumpkins)


#5 Plot the relationship between the estimated weight and actual weight of the pumpkins. 
library(ggplot2)
unique(Pumpkins$weight_class)
ggplot(aes(x = est_weight, y = weight_lbs, color = weight_class), data = Pumpkins) +#Using diamonds data and map carat to x-axis and price to the y-axis
  geom_point(alpha = 0.15, size = 0.01) +
  scale_color_manual(
    values = c(
      light = "#F8766D",   
      medium = "#00BFC4",  
      heavy = "#7CAE00"    
    )
  ) +
  labs(
    x = "Est weight(lbs)",
    y = "Actual weight(lbs)",
    title = "Pumpkin Weight Class Visualization"
  )#done

# 6  Filter the data to include only pumpkins from three countries ofyour choosing and save this filtered data set to your computer incsv format.
library(dplyr)
Pumpkins_n_country <- Pumpkins %>%
  filter(country =="Japan"|
           country=="Italy"|
           country =="Canada")
write.csv(Pumpkins_n_country, "pumpkins_3countries.csv", row.names = FALSE)

#7.1 Summarise your filtered dataset from question 6 ,a. Identify mean weight of pumpkins for each of your three countries — which is highest?
mean_by_country <- Pumpkins_n_country %>%
  group_by(country) %>%
  summarise(mean_weight = mean(weight_lbs, na.rm = TRUE)) %>% # Calculate mean pumpkin weight for each country (ignoring missing values)
  arrange(desc(mean_weight))

mean_by_country

#7.2 b. Identify mean weight for each variety per country — which variety is lowest and where?
lowest_variety_each_country <- Pumpkins_n_country %>%
  group_by(country, variety) %>% ## Group by both country and variety to compute mean per combination
  summarise(mean_weight = mean(weight_lbs, na.rm = TRUE), .groups = "drop") %>% ## Summarise mean weight for each variety within each country
  slice_min(mean_weight, n = 1, by = country) # # Select the variety with the lowest mean weight *within each country*

lowest_variety_each_country
#8 Using filtered data from Q6, draw a boxplot of pumpkin weight distributions (lbs or kg) for the three countries, label axes clearly, save figure.
library(ggplot2)
ggplot(aes(x = country, y= weight_lbs, fill = country), data = Pumpkins_n_country) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "Pumpkin Weight Distribution Across Three Countries",
    x = "Country",
    
    y = "Pumpkin Weight (lbs)"
  ) 

#9 Redraw Q8 using facet plot, with each pumpkin variety split into sub-plots. Save this plot.

library(ggplot2)

ggplot(Pumpkins_n_country, aes(x = country, y = weight_lbs, fill = country)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~ variety, scales = "free_y") +   # 🔥 每个 variety 单独成图
  labs(
    title = "Pumpkin Weight Distribution by Variety and Country",
    x = "Country",
    y = "Pumpkin Weight (lbs)"
  ) 


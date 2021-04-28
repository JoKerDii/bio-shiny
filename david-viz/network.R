library(network)
library(networkD3)
library(dplyr)
library(tidyverse)

mynet <- read.csv("D:/git/bio-shiny/david-viz/data/mybindeddf4networks.csv")
mynet$negLogPvalue <- -log(mynet$PValue)
mynet$negLogPvalue <- (mynet$negLogPvalue - mean(mynet$negLogPvalue))/(max(mynet$negLogPvalue) - min(mynet$negLogPvalue))

From <- mynet %>%
  distinct(Category) %>%
  rename(label = Category)
To <- mynet %>%
  distinct(Term) %>%
  rename(label = Term)

mynodes <- full_join(From, To, by = "label")
mynodes <- mynodes %>% rowid_to_column("id")

myper_route <- mynet %>%  
  select(Category, Term, negLogPvalue) 

myedges <- myper_route %>% 
  left_join(mynodes, by = c("Category" = "label")) %>% 
  rename(from = id)

myedges <- myedges %>% 
  left_join(mynodes, by = c("Term" = "label")) %>% 
  rename(to = id)

myedges <- select(myedges, from, to, negLogPvalue)
# myroutes_network <- network(myedges, vertex.attr = mynodes, matrix.type = "edgelist", ignore.eval = FALSE)

nodes_d3 <- mutate(mynodes, id = id - 1)
edges_d3 <- mutate(myedges, from = from - 1, to = to - 1)


forceNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
             NodeID = "label", Group = "id", Value = "negLogPvalue", 
             opacity = 1, fontSize = 16, zoom = FALSE)

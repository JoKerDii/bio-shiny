library(dplyr)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(highcharter)
library(stringr)
library(httr)
library(rdrop2)
library(readr)

local.path.to.tables = "data/"
chart <- read_tsv(paste(local.path.to.tables,"chart.txt", sep = ""))
table <- read_tsv(paste(local.path.to.tables,"table.txt", sep = ""))
# > unique(chart$Category)
# [1] "UP_KEYWORDS"      "GOTERM_MF_DIRECT" "GOTERM_CC_DIRECT"
# [4] "UP_SEQ_FEATURE"   "GOTERM_BP_DIRECT" "INTERPRO"        
# [7] "SMART"            "KEGG_PATHWAY"     "BIOCARTA"        
# [10] "PIR_SUPERFAMILY"  "OMIM_DISEASE"     "BBID"            
# [13] "COG_ONTOLOGY" 


# strsplit(table$GOTERM_BP_DIRECT, ",")

# gsub("^.*?~", "", table_tmp$GOTERM_BP_DIRECT_S)

# table %>%
#   rename(Gene_Name = `Gene Name`) %>%
#   separate_rows(sep = ",", GOTERM_BP_DIRECT, convert = TRUE) %>%
#   mutate(GOTERM_BP_DIRECT = as.factor(trimws(gsub("^.*?~", "",GOTERM_BP_DIRECT)))) %>%
#   filter(GOTERM_BP_DIRECT != "")

BIOCARTA <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(BIOCARTA != "NA") %>%
  select(ID, Gene_Name, BIOCARTA) %>%
  mutate(Group = "BIOCARTA") %>%
  rename(Terms = BIOCARTA) %>%
  as.data.frame()


COG_ONTOLOGY <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(COG_ONTOLOGY != "NA") %>%
  select(ID, Gene_Name, COG_ONTOLOGY) %>%
  filter(COG_ONTOLOGY != "") %>%
  mutate(Group = "COG_ONTOLOGY") %>%
  rename(Terms = COG_ONTOLOGY) %>%
  as.data.frame()


GOTERM_BP_DIRECT <-table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(GOTERM_BP_DIRECT != "NA") %>%
  select(ID, Gene_Name, GOTERM_BP_DIRECT)  %>%
  separate_rows(sep = ",", GOTERM_BP_DIRECT, convert = TRUE) %>%
  mutate(GOTERM_BP_DIRECT = as.factor(gsub("^.*?~", "",GOTERM_BP_DIRECT))) %>%
  filter(GOTERM_BP_DIRECT != "") %>%
  mutate(Group = "GOTERM_BP_DIRECT") %>%
  rename(Terms = GOTERM_BP_DIRECT) %>%
  as.data.frame()

GOTERM_CC_DIRECT <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(GOTERM_CC_DIRECT != "NA") %>%
  select(ID, Gene_Name, GOTERM_CC_DIRECT)  %>%
  separate_rows(sep = ",", GOTERM_CC_DIRECT, convert = TRUE) %>%
  mutate(GOTERM_CC_DIRECT = as.factor(gsub("^.*?~", "",GOTERM_CC_DIRECT))) %>%
  filter(GOTERM_CC_DIRECT != "")%>%
  mutate(Group = "GOTERM_CC_DIRECT") %>%
  rename(Terms = GOTERM_CC_DIRECT) %>%
  as.data.frame()

GOTERM_MF_DIRECT <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(GOTERM_MF_DIRECT != "NA") %>%
  select(ID, Gene_Name, GOTERM_MF_DIRECT)  %>%
  separate_rows(sep = ",", GOTERM_MF_DIRECT, convert = TRUE) %>%
  mutate(GOTERM_MF_DIRECT = as.factor(gsub("^.*?~", "",GOTERM_MF_DIRECT))) %>%
  filter(GOTERM_MF_DIRECT != "")%>%
  mutate(Group = "GOTERM_MF_DIRECT") %>%
  rename(Terms = GOTERM_MF_DIRECT) %>%
  as.data.frame()


INTERPRO <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(INTERPRO != "NA") %>%
  select(ID, Gene_Name, INTERPRO)  %>%
  separate_rows(sep = ",IPR", INTERPRO, convert = TRUE) %>%
  mutate(INTERPRO = as.factor(gsub("^.*?:", "",INTERPRO))) %>%
  filter(INTERPRO != "")%>%
  mutate(Group = "INTERPRO") %>%
  rename(Terms = INTERPRO) %>%
  as.data.frame()


KEGG_PATHWAY <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(KEGG_PATHWAY != "NA") %>%
  select(ID, Gene_Name, KEGG_PATHWAY) %>%
  separate_rows(sep = ",hsa", KEGG_PATHWAY, convert = TRUE) %>%
  mutate(KEGG_PATHWAY = as.factor(gsub("^.*?:", "",KEGG_PATHWAY))) %>%
  filter(KEGG_PATHWAY != "")%>%
  mutate(Group = "KEGG_PATHWAY") %>%
  rename(Terms = KEGG_PATHWAY) %>%
  as.data.frame()


OMIM_DISEASE <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(OMIM_DISEASE != "NA") %>%
  select(ID, Gene_Name, OMIM_DISEASE) %>%
  separate_rows(sep = ",.*?~", OMIM_DISEASE, convert = TRUE) %>%
  mutate(OMIM_DISEASE = as.factor(gsub("^.*?~", "",OMIM_DISEASE))) %>%
  filter(OMIM_DISEASE != "")%>%
  mutate(Group = "OMIM_DISEASE") %>%
  rename(Terms = OMIM_DISEASE) %>%
  as.data.frame()

PIR_SUPERFAMILY <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(PIR_SUPERFAMILY != "NA") %>%
  select(ID, Gene_Name, PIR_SUPERFAMILY) %>%
  separate_rows(sep = ",P.*?:", PIR_SUPERFAMILY, convert = TRUE) %>%
  mutate(PIR_SUPERFAMILY = as.factor(gsub("^.*?:", "",PIR_SUPERFAMILY))) %>%
  filter(PIR_SUPERFAMILY != "")%>%
  mutate(Group = "PIR_SUPERFAMILY") %>%
  rename(Terms = PIR_SUPERFAMILY) %>%
  as.data.frame()


SMART <- table %>%
  rename(Gene_Name = `Gene Name`) %>%
  filter(SMART != "NA") %>%
  select(ID, Gene_Name, SMART) %>%
  separate_rows(sep = ",S.*?:", SMART, convert = TRUE) %>%
  mutate(SMART = as.factor(gsub("^.*?:", "",SMART))) %>%
  filter(SMART != "")%>%
  mutate(Group = "SMART") %>%
  rename(Terms = SMART) %>%
  as.data.frame()


mybindeddf <- bind_rows(BIOCARTA, COG_ONTOLOGY, GOTERM_BP_DIRECT, GOTERM_CC_DIRECT,
                        GOTERM_MF_DIRECT, INTERPRO, KEGG_PATHWAY, OMIM_DISEASE, 
                        PIR_SUPERFAMILY, SMART)

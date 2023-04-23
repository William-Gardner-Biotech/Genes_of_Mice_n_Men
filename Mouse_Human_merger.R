setwd("~/Desktop")

library(tidyverse)
library(dplyr)

Mus_v_humo <- read_csv(file = "Mouse_Human_clean.csv")

#view(Mus_v_humo)

mice <- read_tsv(file = "mouse_brain_category_rna_cerebral.tsv")

mouse_rough <- mice

#view(mouse_rough)

humans <- read_tsv(file = "rna_tissue_gtex.tsv")

humans = humans[-c(4,5)]

human_cereb <- filter(humans, Tissue == "cerebellum")

colnames(human_cereb)[4] = "Human_nTPM"

human_cereb <- human_cereb[-c(3)]

colnames(human_cereb)[1:2] = c("Human_ENSG", "Gene_name")

#view(human_cereb)

#view(human_cereb)

#view(Mus_v_humo)

blast_Mus_v_humo <- unique(Mus_v_humo)

colnames(blast_Mus_v_humo) = c("Mouse_gene", "Human_gene")

#view(blast_Mus_v_humo)

mice <- mice[c(1:17,93)]

mice_awk <- mice[-c(2,6,7,12:16)]

colnames(mice_awk)[5] = "Protein_desc"

view(mice_awk)

merge_1 <- left_join(blast_Mus_v_humo, human_cereb, by=c("Human_gene"= "Gene_name"))

merge_1 %>% drop_na("Human_ENSG") 

merge_1 <- filter(merge_1, Human_nTPM > 1)

view(merge_1)

merge_2 <- left_join(merge_1, mice_awk, by=c("Mouse_gene"= "Gene"))

merge_2 <- merge_2 %>% relocate(13, .after = 4)

merge_2 <- merge_2 %>% drop_na("Protein_desc") 

colnames(merge_2)[5] = "Mouse_nTPM"

colnames(merge_2)[13] = "RNA_tissue_specificity"

merge_2 <- merge_2[-c(3,6,8)]

view(merge_2)

write_csv( merge_2, file = "Final_table.csv")

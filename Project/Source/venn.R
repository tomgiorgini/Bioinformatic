# 1) librerie
library(dplyr)
library(stringr)
setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic/Project")

C_AD = read.table('Results/C_AD/DEG.txt', header = T, sep = "\t", quote = "", check.names=F)
C_MCI = read.table('Results/C_MCI/DEG.txt', header = T, sep = "\t", quote = "", check.names=F)
AD_MCI = read.table('Results/AD_MCI/DEG.txt', header = T, sep = "\t", quote = "", check.names=F)

# 2) funzione per parsare i file venn_*.txt
parse_venn <- function(file) {
  lines <- readLines(file)
  df <- tibble(GeneSymbol = character(), set = character())
  for(line in lines) {
    # split in "chiave" e "lista geni"
    parts    <- str_split_fixed(line, ":", n = 2)
    set_info <- str_remove_all(parts[1], "\\[|\\]")      # es. "C_AD and C_MCI"
    genes    <- str_split(parts[2], ",")[[1]] %>% str_trim()
    df <- bind_rows(df, tibble(GeneSymbol = genes, set = set_info))
  }
  df
}

# 3) leggi e parsare UP e DOWN
venn_up   <- parse_venn("venn_up.txt")   %>% mutate(direction = "UP")
venn_down <- parse_venn("venn_down.txt") %>% mutate(direction = "DOWN")

# unisco
venn_info <- bind_rows(venn_up, venn_down)


# 4) metti le tue tre matrici in una lista  
#    (assumendo che C_AD, C_MCI e AD_MCI siano già in memoria come data.frame)
lst <- list(C_AD = C_AD, C_MCI = C_MCI, AD_MCI = AD_MCI)

# aggiungo colonna "study" e le concateno
all_df <- bind_rows(lst, .id = "study")


# 5) unisco i dati originali con l’informazione di 'set'  
#    (faccio un inner join su GeneSymbol + direction,
#     in modo da mantenere solo i geni che hai passato nei venn_*.txt)


final <- all_df %>%
  inner_join(venn_info, by = c("GeneSymbol", "direction")) %>%
  select(GeneSymbol, set, direction) %>%
  distinct()

final <- final %>% arrange(set,direction)
# 6) (opzionale) esporta
write.table(final, "genes_sets.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(final$GeneSymbol[final$set == 'C_AD' & final$direction == 'UP'], "C_AD_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_AD' & final$direction == 'DOWN'], "C_AD_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_MCI' & final$direction == 'UP'], "C_MCI_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_MCI' & final$direction == 'DOWN'], "C_MCI_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'AD_MCI' & final$direction == 'UP'], "AD_MCI_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'AD_MCI' & final$direction == 'DOWN'], "AD_MCI_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_AD and C_MCI' & final$direction == 'UP'], "C_AD_and_C_MCI_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_AD and C_MCI' & final$direction == 'DOWN'], "C_AD_and_C_MCI_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = F)
write.table(final$GeneSymbol[final$set == 'C_MCI and AD_MCI' & final$direction == 'UP'], "C_MCI_and_AD_MCI_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_MCI and AD_MCI' & final$direction == 'DOWN'], "C_MCI_and_AD_MCI_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_AD and C_MCI and AD_MCI' & final$direction == 'UP'], "C_AD_and_C_MCI_and_AD_MCI_UP.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)
write.table(final$GeneSymbol[final$set == 'C_AD and C_MCI and AD_MCI' & final$direction == 'DOWN'], "C_AD_and_C_MCI_and_AD_MCI_DOWN.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,col.names = F)

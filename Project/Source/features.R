# 1) librerie
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
setwd("C:/Users/tomma/OneDrive/Desktop/Bioinformatic/Project")


# 2) leggi i metadata
meta <- read_tsv("Data/metadata.txt", col_names = TRUE) %>%
  # 1) rinomina le colonne “with :ch1” in nomi semplici
  rename(
    age       = `age:ch1`,
    gender    = `gender:ch1`,
    ethnicity = `ethnicity:ch1`,
    sample_status    = `status`,
    disease_status = 'status:ch1'
  ) %>%
  # 2) ora mutate sui nomi nuovi esiste davvero
  mutate(
    age       = as.numeric(age),
    gender    = as.factor(gender),
    ethnicity = as.factor(ethnicity),
    status    = as.factor(sample_status)
  )

meta <- meta %>%
  mutate(
    age_group = cut(age,
                    breaks = seq(55, 80, by = 5),
                    right  = FALSE,
                    labels = paste(seq(55,75,by=5), seq(59,79,by=5), sep = "-"))
  )

# 3) crea la cartella di output
out_dir <- "plots"
if (!dir.exists(out_dir)) dir.create(out_dir)

# 4) Histogram per age con blocchi di 5 anni
p_age <- ggplot(meta, aes(x = age)) +
  geom_histogram(binwidth = 5, color = "black", fill = "steelblue", alpha = 0.7) +
  scale_x_continuous(breaks = seq(floor(min(meta$age, na.rm=TRUE)),
                                  ceiling(max(meta$age, na.rm=TRUE)),
                                  by = 5)) +
  labs(title = "Distribuzione dell'età (blocchi di 5 anni)",
       x = "Età (anni)",
       y = "Frequenza") +
  theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = file.path(out_dir, "age_histogram_5yr.png"),
  plot     = p_age,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# 5) funzione helper per pie chart con conteggi (allineamento corretto delle etichette)
make_pie_count <- function(df, var, filename, title) {
  counts <- df %>%
    count(!!sym(var)) %>%
    arrange(desc(n))
  
  p <- ggplot(counts, aes(x = "", y = n, fill = !!sym(var))) +
    geom_col(width = 1, color = "white", show.legend = TRUE) +
    coord_polar(theta = "y") +
    geom_text(aes(label = n), position = position_stack(vjust = 0.7), size = 3) +
    labs(title = title, fill = NULL) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "right",
      plot.title       = element_text(hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(out_dir, filename),
    plot     = p,
    width    = 5,
    height   = 5,
    dpi      = 300
  )
}

# 6) genera e salva i pie chart con conteggi
make_pie_count(meta, "gender",         "gender_pie_counts.png",         "Genere")
make_pie_count(meta, "ethnicity",      "ethnicity_pie_counts.png",      "Etnia")
make_pie_count(meta, "disease_status", "disease_status_pie_counts.png", "Disease Status")
library(tidyverse)
library(readxl)
library(cowplot)
library(MESS)
library(rstatix)
library(ggpubr)

# ----------------------------------------------------------------------
# 1. Configuración y Tema Nature
# ----------------------------------------------------------------------
out_dir <- "C:/Users/santi/OneDrive/Escritorio/BVDV Imagenes/Manuscrito/Signos clinicos"
figures_dir <- file.path(out_dir, "figuras")
if(!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

palette_nature <- c(
  "Control"    = "#4D4D4D",
  "25 µg"      = "#004B87",
  "50 µg"      = "#E69F00",
  "100 µg"     = "#2ECC71",
  "Commercial" = "#D55E00"
)

shapes_nature <- c(
  "Control"    = 16,
  "25 µg"      = 15,
  "50 µg"      = 17,
  "100 µg"     = 25,
  "Commercial" = 18
)

group_levels <- c("Control", "25 µg", "50 µg", "100 µg", "Commercial")

theme_nature <- function(base_size = 11) {
  theme_cowplot(font_size = base_size, font_family = "arial") +
    theme(
      plot.title           = element_text(face = "plain", size = base_size + 1, hjust = 0.5),
      axis.title           = element_text(size = base_size),
      axis.text            = element_text(size = base_size - 1, color = "grey15"),
      legend.title         = element_blank(),
      legend.position      = "top",
      legend.justification = "center",
      legend.text          = element_text(size = base_size - 1),
      panel.grid.major.y   = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank()
    )
}

# ----------------------------------------------------------------------
# 2. Carga y Preparación de Datos
# ----------------------------------------------------------------------
excel_path <- file.path(out_dir, "Diarrea.xlsx")
df <- read_excel(excel_path)
colnames(df) <- str_trim(colnames(df))

df_clean <- df %>%
  mutate(
    Score_raw = as.numeric(gsub(",", ".", as.character(Score))),
    Score = case_when(
      Score_raw == 1.2 ~ 1.5,
      Score_raw == 0.1 ~ 0.5,
      TRUE             ~ Score_raw
    ),
    Tiempo = as.numeric(Tiempo),
    DIIO = factor(DIIO),
    Vacuna_Clean = str_trim(Vacuna),
    Grupo_Factor = case_when(
      Vacuna_Clean == "Control (-)"       ~ "Control",
      Vacuna_Clean == "rE2 25ug"          ~ "25 µg",
      Vacuna_Clean == "rE2 50ug"          ~ "50 µg",
      Vacuna_Clean == "rE2 100ug"         ~ "100 µg",
      Vacuna_Clean == "Cattle master C+"  ~ "Commercial",
      TRUE                                ~ Vacuna_Clean
    ),
    Grupo_Factor = factor(Grupo_Factor, levels = group_levels)
  )

# Resumen Longitudinal
df_summary <- df_clean %>%
  group_by(Tiempo, Grupo_Factor) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    se_score   = sd(Score, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  )

# Cálculo de AUC por animal
df_auc <- df_clean %>%
  group_by(DIIO, Grupo_Factor) %>%
  arrange(Tiempo) %>%
  summarise(
    AUC = MESS::auc(Tiempo, Score),
    .groups = "drop"
  )

# ----------------------------------------------------------------------
# 3. ANÁLISIS ESTADÍSTICO Y SALIDA DE RESULTADOS A CONSOLA
# ----------------------------------------------------------------------

cat("\n======================================================================\n")
cat(" 1. ESTADÍSTICA DESCRIPTIVA DE AUC POR GRUPO\n")
cat("======================================================================\n")
auc_descriptive <- df_auc %>%
  group_by(Grupo_Factor) %>%
  summarise(
    N       = n(),
    Median  = median(AUC, na.rm = TRUE),
    IQR     = IQR(AUC, na.rm = TRUE),
    Mean    = mean(AUC, na.rm = TRUE),
    SD      = sd(AUC, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(auc_descriptive))

cat("\n======================================================================\n")
cat(" 2. PRUEBA GLOBAL DE KRUSKAL-WALLIS PARA AUC\n")
cat("======================================================================\n")
kw_res <- df_auc %>% kruskal_test(AUC ~ Grupo_Factor)
print(as.data.frame(kw_res))

cat("\n======================================================================\n")
cat(" 3. POST-HOC: PRUEBA DE DUNN (TODAS LAS COMPARACIONES PAREADAS)\n")
cat("======================================================================\n")
dunn_all <- df_auc %>%
  dunn_test(AUC ~ Grupo_Factor, p.adjust.method = "fdr")

# Mostrar todas las comparaciones en consola
print(as.data.frame(dunn_all))

cat("\n--- Comparaciones significativas (p.adj < 0.05) ---\n")
dunn_sig <- dunn_all %>% filter(p.adj < 0.05)
if(nrow(dunn_sig) > 0) {
  print(as.data.frame(dunn_sig))
} else {
  cat("No se encontraron diferencias estadísticamente significativas (p.adj >= 0.05).\n")
}
cat("======================================================================\n\n")

# Preparar dataset filtrado para la figura (solo p.adj < 0.05)
if(nrow(dunn_sig) > 0) {
  stat_test_plot <- dunn_sig %>%
    add_y_position(step.increase = 0.12) %>%
    mutate(p.adj.signif = ifelse(p.adj < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p.adj))))
} else {
  stat_test_plot <- dunn_sig
}

# ----------------------------------------------------------------------
# 4. Gráficos
# ----------------------------------------------------------------------
p_ribbon <- ggplot(df_summary, aes(x = Tiempo, y = mean_score, color = Grupo_Factor, fill = Grupo_Factor)) +
  geom_ribbon(aes(ymin = pmax(0, mean_score - se_score), ymax = mean_score + se_score), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Grupo_Factor), size = 2.2) +
  scale_color_manual(values = palette_nature) +
  scale_fill_manual(values = palette_nature) +
  scale_shape_manual(values = shapes_nature) +
  scale_y_continuous(
    limits = c(0, 3.1), 
    breaks = 0:3, 
    labels = c("0: Normal", "1: Mild", "2: Watery", "3: Bloody")
  ) +
  scale_x_continuous(breaks = seq(0, max(df_clean$Tiempo, na.rm = TRUE), by = 2)) +
  labs(
    x = "Days after challenge", 
    y = "Diarrhea Score (Mean ± SE)", 
    title = "Clinical Diarrhea Course"
  ) +
  theme_nature()

p_auc <- ggplot(df_auc, aes(x = Grupo_Factor, y = AUC, fill = Grupo_Factor, color = Grupo_Factor)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.5, linewidth = 0.6) +
  geom_jitter(aes(shape = Grupo_Factor), width = 0.15, size = 2, alpha = 0.8) +
  scale_color_manual(values = palette_nature) +
  scale_fill_manual(values = palette_nature) +
  scale_shape_manual(values = shapes_nature) +
  labs(
    x = "Vaccination Group", 
    y = "Cumulative Diarrhea Score (AUC)", 
    title = "Total Disease Burden"
  ) +
  theme_nature() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

if(nrow(stat_test_plot) > 0) {
  p_auc <- p_auc + 
    stat_pvalue_manual(
      stat_test_plot,
      label = "p.adj.signif",
      tip.length = 0.015,
      size = 3.3,
      color = "grey20"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
}

# ----------------------------------------------------------------------
# 5. Combinación de Paneles y Guardado
# ----------------------------------------------------------------------
p_combined <- plot_grid(
  p_ribbon, p_auc,
  labels = c("A", "B"),
  label_size = 14,
  rel_widths = c(1.3, 1)
)

print(p_combined)

ggsave(file.path(figures_dir, "figure_diarrhea_ribbon.tiff"), plot = p_ribbon, width = 5, height = 3.5, units = "in", dpi = 600, compression = "lzw")
ggsave(file.path(figures_dir, "figure_diarrhea_auc_boxplot.tiff"), plot = p_auc, width = 2.1, height = 3, units = "in", dpi = 600, compression = "lzw")
ggsave(file.path(figures_dir, "figure_diarrhea_combined.tiff"), plot = p_combined, width = 10.5, height = 4.8, units = "in", dpi = 600, compression = "lzw")
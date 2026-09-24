# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# ANÁLISIS SERONEUTRALIZACIÓN (log2) - BOXPLOT PANEL 2x2 & EXCEL
#
# Autor: Santiago Salazar
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Cargar y verificar paquetes
# ------------------------------------------------------------
required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "openxlsx",
  "lme4",
  "lmerTest",
  "emmeans",
  "patchwork",
  "cowplot",
  "purrr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(openxlsx)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(patchwork)
  library(cowplot)
  library(purrr)
})

set.seed(20260620)


# ------------------------------------------------------------
# 2. Localización del Proyecto y Directorios
# ------------------------------------------------------------
find_project_root <- function() {
  current_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (current_dir != dirname(current_dir)) {
    rproj_files <- list.files(current_dir, pattern = "\\.Rproj$", full.names = TRUE)
    if (length(rproj_files) > 0) return(current_dir)
    current_dir <- dirname(current_dir)
  }
  return(getwd())
}

project_dir <- find_project_root()

data_dir    <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
tables_dir  <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
models_dir  <- file.path(results_dir, "models")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Lectura y Preparación de Datos
# ------------------------------------------------------------
excel_files <- list.files(data_dir, pattern = "^Neutralization.*\\.(xlsx|xls)$", full.names = TRUE)
excel_path <- if (length(excel_files) > 0) excel_files[1] else file.path(data_dir, "Neutralization.xlsx")

if (!file.exists(excel_path)) {
  stop(paste0("No se encontró el archivo de datos en: ", excel_path))
}

# Modificado para reflejar las etiquetas de tratamiento (5 niveles)
grupo_levels <- c("Control", "25 µg", "50 µg", "100 µg", "Comm.")

datos_titulos <- read_excel(excel_path) %>%
  mutate(
    Grupo = trimws(Grupo),
    # Mapeo flexible de nombres si vienen como '25ug', 'Comercial', etc.
    Grupo = case_when(
      Grupo %in% c("Control") ~ "Control",
      Grupo %in% c("25ug", "25 µg rE2", "25ug rE2") ~ "25 µg",
      Grupo %in% c("50ug", "50 µg rE2", "50ug rE2") ~ "50 µg",
      Grupo %in% c("100ug", "100 µg rE2", "100ug rE2") ~ "100 µg",
      Grupo %in% c("Comercial", "Comm.", "Commercial") ~ "Comm.",
      TRUE ~ Grupo
    ),
    Grupo = factor(Grupo, levels = grupo_levels),
    Virus = factor(trimws(Virus)),
    Animal = factor(Animal),
    Log2_Titulo = log2(Titulo)
  )

# Subset para inferencia estadística (sin Control)
datos_stats <- datos_titulos %>%
  filter(Grupo != "Control") %>%
  mutate(Grupo = factor(Grupo, levels = c("25 µg", "50 µg", "100 µg", "Comm.")))


# ------------------------------------------------------------
# 4. Estilos y Mapeo de Virus
# ------------------------------------------------------------
palette_nature <- c(
  "Control"    = "#595959",
  "25 µg"  = "#004B87",
  "50 µg"  = "#E69F00",
  "100 µg" = "#2ECC71",
  "Comm."      = "#D55E00"
)

shapes_nature <- c(
  "Control"    = 21,
  "25 µg"  = 22,
  "50 µg"  = 24,
  "100 µg" = 25,
  "Comm."      = 23
)

theme_nature <- function(base_size = 10) {
  theme_cowplot(font_size = base_size, font_family = "arial") +
    theme(
      plot.title         = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      axis.title         = element_text(size = base_size, face = "plain"),
      axis.title.y       = element_text(margin = margin(r = 6)),
      axis.title.x       = element_text(margin = margin(t = 6)),
      axis.text          = element_text(size = base_size - 1, color = "black"),
      axis.line          = element_line(color = "black", linewidth = 0.5),
      legend.title       = element_blank(),
      legend.position    = "bottom",
      legend.text        = element_text(size = base_size - 1),
      legend.background  = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

virus_map <- list(
  "B" = list(name = "1b", title = "BVDV-1b"),
  "C" = list(name = "1c", title = "BVDV-1c"),
  "D" = list(name = "1d", title = "BVDV-1d"),
  "E" = list(name = "1e", title = "BVDV-1e")
)


# ------------------------------------------------------------
# 5. Ajuste del Modelo Lineal Mixto (LMM)
# ------------------------------------------------------------

lmm_model <- lmer(Log2_Titulo ~ Grupo * Virus + (1 | Animal), data = datos_stats)
saveRDS(lmm_model, file.path(models_dir, "Neutralization_LMM.rds"))

anova_df <- as.data.frame(anova(lmm_model))
anova_df <- cbind(Efecto = rownames(anova_df), anova_df)

emms_res <- emmeans(lmm_model, ~ Grupo | Virus)
emms_df  <- as.data.frame(emms_res)

pairs_res <- pairs(emms_res, adjust = "BH")
pairs_df  <- as.data.frame(pairs_res)


# ------------------------------------------------------------
# 6. Función de Generación de Boxplots por Subgenotipo
# ------------------------------------------------------------
generate_subgenotype_boxplot <- function(cepa_key) {
  meta <- virus_map[[cepa_key]]
  
  plot_data <- datos_titulos %>%
    filter(Virus == cepa_key) %>%
    filter(!is.na(Log2_Titulo))
  
  stats_virus <- pairs_df %>%
    filter(Virus == cepa_key)
  
  max_y <- max(plot_data$Log2_Titulo, na.rm = TRUE)
  min_y <- min(plot_data$Log2_Titulo, na.rm = TRUE)
  range_y <- max_y - min_y
  if (range_y == 0) range_y <- 1
  
  # Alinear búsquedas con las posiciones de x (1=Control, 2=25, 3=50, 4=100, 5=Comm.)
  comp_metadata <- list(
    list(g1 = "25 µg",  x1 = 2, x2 = 5, row = 1, pattern = "25 µg - Comm.|Comm. - 25 µg"),
    list(g1 = "50 µg",  x1 = 3, x2 = 5, row = 2, pattern = "50 µg - Comm.|Comm. - 50 µg"),
    list(g1 = "100 µg", x1 = 4, x2 = 5, row = 3, pattern = "100 µg - Comm.|Comm. - 100 µg")
  )
  
  lines_and_labels <- map_df(comp_metadata, function(m) {
    p_val <- stats_virus %>%
      filter(grepl(m$pattern, contrast)) %>%
      pull(p.value)
    
    if (length(p_val) == 0) p_val <- 1.0
    
    p_formatted <- if (p_val < 0.001) "p < 0.001" else if (p_val >= 0.05) "ns" else paste0("p = ", sprintf("%.3f", p_val))
    y_bar <- max_y + (range_y * 0.10 * m$row)
    
    data.frame(
      x = m$x1, xend = m$x2, y = y_bar,
      x_text = (m$x1 + m$x2) / 2, y_text = y_bar + (range_y * 0.03),
      p_label = p_formatted
    )
  })
  
  ylim_max <- max(lines_and_labels$y_text) + (range_y * 0.05)
  
  p <- ggplot(plot_data, aes(x = Grupo, y = Log2_Titulo, color = Grupo, fill = Grupo)) +
    geom_boxplot(
      width = 0.5, 
      outlier.shape = NA, 
      alpha = 0.15, 
      linewidth = 0.4
    ) +
    geom_point(
      aes(shape = Grupo),
      position = position_jitterdodge(
        jitter.width = 0.35, 
        jitter.height = 0, 
        dodge.width = 0.5,
        seed = 20260620
      ),
      size = 1.6, 
      alpha = 0.85,
      stroke = 0.3
    ) +
    geom_segment(data = lines_and_labels, aes(x = x, xend = xend, y = y, yend = y), 
                 color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_segment(data = lines_and_labels, aes(x = x, xend = x, y = y, yend = y - (range_y * 0.02)), 
                 color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_segment(data = lines_and_labels, aes(x = xend, xend = xend, y = y, yend = y - (range_y * 0.02)), 
                 color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_text(data = lines_and_labels, aes(x = x_text, y = y_text, label = p_label), 
              color = "black", size = 2.2, vjust = 0, fontface = "plain", inherit.aes = FALSE) +
    scale_color_manual(values = palette_nature, drop = FALSE) +
    scale_fill_manual(values = palette_nature, drop = FALSE) +
    scale_shape_manual(values = shapes_nature, drop = FALSE) +
    coord_cartesian(ylim = c(0, 18)) +
    scale_y_continuous(breaks = seq(0, 16, by = 2))+
    labs(
      title = meta$title,
      x = NULL, 
      y = expression("Neutralizing Antibody Titer (" * log[2] * ")")
    ) +
    theme_nature(10) +
    theme(
      plot.margin = margin(t = 5, r = 8, b = 5, l = 5)
    )
  
  return(p)
}


# ------------------------------------------------------------
# 7. Integrar y Guardar Panel 2x2
# ------------------------------------------------------------
cepas <- names(virus_map)
list_boxplots <- map(cepas, generate_subgenotype_boxplot)

# Ensamblado del panel 2x2 con patchwork
panel_2x2_boxplot <- wrap_plots(list_boxplots, ncol = 2, nrow = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_path <- file.path(figures_dir, "Seroneutralization_BVDV")

ggsave(paste0(fig_path, ".png"), panel_2x2_boxplot, width = 6, height = 7, dpi = 300)
ggsave(paste0(fig_path, ".pdf"), panel_2x2_boxplot, width = 6, height = 7)
ggsave(paste0(fig_path, ".tiff"), panel_2x2_boxplot, width = 6, height = 7, dpi = 600, compression = "lzw")

message("Figura Boxplot 2x2 generada en formatos PNG, PDF y TIFF.")


# ------------------------------------------------------------
# 8. Exportación de Resultados a Un Solo Archivo Excel
# ------------------------------------------------------------
output_excel <- file.path(tables_dir, "Seroneutralization_Analysis.xlsx")
wb <- createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  if (is.null(data) || nrow(data) == 0) {
    writeData(wb, sheet = sheet_name, x = data.frame(Information = "No data available"))
  } else {
    writeData(wb, sheet = sheet_name, x = data)
  }
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
}

# Resumen Descriptivo
resumen_stats <- datos_titulos %>%
  group_by(Virus, Grupo) %>%
  summarise(
    N = n(),
    Media_Log2 = mean(Log2_Titulo, na.rm = TRUE),
    SD_Log2 = sd(Log2_Titulo, na.rm = TRUE),
    SEM_Log2 = SD_Log2 / sqrt(N),
    Mediana_Log2 = median(Log2_Titulo, na.rm = TRUE),
    IQR_Log2 = IQR(Log2_Titulo, na.rm = TRUE),
    .groups = "drop"
  )

readme_table <- data.frame(
  Section = c("Input Data", "Descriptive Stats", "LMM ANOVA", "EMMeans", "Pairwise Contrasts"),
  Description = c(
    "Log2 neutralizing antibody titers dataset",
    "Summary stats by virus and group (Mean, SD, SEM, Median, IQR)",
    "Linear Mixed-Effects Model Type III ANOVA table",
    "Estimated Marginal Means for treatment groups",
    "Pairwise comparisons (BH adjusted p-values)"
  )
)

add_sheet(wb, "README", readme_table)
add_sheet(wb, "Datos_Procesados", datos_titulos)
add_sheet(wb, "Resumen_Estadistico", resumen_stats)
add_sheet(wb, "ANOVA_LMM", anova_df)
add_sheet(wb, "EMMeans", emms_df)
add_sheet(wb, "Comparaciones_Pares", pairs_df)

saveWorkbook(wb, output_excel, overwrite = TRUE)
message("Análisis exportado exitosamente a: ", output_excel)
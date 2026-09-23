# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# ANÁLISIS NO PARAMÉTRICO DE TÍTULOS DE ANTICUERPOS (log2)
#
# Autor: Santiago Salazar
# ============================================================

options(stringsAsFactors = FALSE)

# 1. Cargar y verificar paquetes
required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "openxlsx",
  "dunn.test",
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
  library(dunn.test)
  library(patchwork)
  library(cowplot)
  library(purrr)
})

set.seed(20260620)


# ============================================================
# 2. LOCALIZACIÓN DEL PROYECTO Y DIRECTORIOS (GITHUB / RPROJ)
# ============================================================

find_project_root <- function() {
  current_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  while (current_dir != dirname(current_dir)) {
    rproj_files <- list.files(current_dir, pattern = "\\.Rproj$", full.names = TRUE)
    if (length(rproj_files) > 0) return(current_dir)
    current_dir <- dirname(current_dir)
  }
  stop("No se encontró el archivo .Rproj del proyecto.")
}

project_dir <- find_project_root()

data_dir    <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
tables_dir  <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 3. LECTURA Y PREPARACIÓN DE DATOS
# ============================================================

excel_path <- file.path(data_dir, "Titulo.xlsx")

if (!file.exists(excel_path)) {
  stop(paste0("No se encontró el archivo:\n", excel_path))
}

datos_titulos <- read_excel(excel_path) %>%
  mutate(
    Grupo = factor(Grupo, levels = c("Control", "25ug", "50ug", "100ug", "Comercial")),
    Animal = factor(Animal),
    Tiempo_Factor = factor(paste0("Día ", Tiempo), levels = c("Día 28", "Día 42")),
    Tiempo_Num = as.numeric(Tiempo),
    Log2_Titulo = log2(as.numeric(Título))
  )


# ============================================================
# 4. PALETA DE COLOR Y TEMA NATURE
# ============================================================

palette_nature <- c(
  "Control"   = "#4D4D4D", # Gris oscuro
  "25ug"      = "#004B87", # Azul oscuro
  "50ug"      = "#E69F00", # Naranja
  "100ug"     = "#2ECC71", # Verde
  "Comercial" = "#D55E00"  # Rojo bermellón
)

theme_nature <- function(base_size = 11) {
  cowplot::theme_cowplot(font_size = base_size, font_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size, hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, color = "grey25"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, color = "grey15"),
      legend.title = element_blank(),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

save_nature <- function(plot, filename, width = 7, height = 3.8) {
  ggsave(file.path(figures_dir, paste0(filename, ".tiff")), plot, width = width, height = height, units = "in", dpi = 600, compression = "lzw")
  ggsave(file.path(figures_dir, paste0(filename, ".pdf")),  plot, width = width, height = height, units = "in")
  ggsave(file.path(figures_dir, paste0(filename, ".png")),  plot, width = width, height = height, units = "in", dpi = 300)
}


# ============================================================
# 5. FIGURA 1: COMPARACIÓN ENTRE GRUPOS EN CADA DÍA (KRUSKAL-WALLIS + DUNN BH)
# ============================================================

metric_metadata <- list(
  "Dia_28" = list(
    filter_expr = quote(Tiempo_Num == 28),
    title = "Endpoint Titer\n(Day 28)",
    ylabel = expression("Endpoint Antibody Titer (" * log[2] * ")")
  ),
  "Dia_42" = list(
    filter_expr = quote(Tiempo_Num == 42),
    title = "Endpoint Titer\n(Day 42)",
    ylabel = expression("Endpoint Antibody Titer (" * log[2] * ")")
  )
)

dunn_results_list <- list()

generate_endpoint_titer_boxplot <- function(meta, day_name) {
  plot_data <- datos_titulos %>% 
    filter(eval(meta$filter_expr)) %>% 
    select(Grupo, Value = Log2_Titulo) %>% 
    filter(!is.na(Value))
  
  stats_data <- plot_data %>% filter(Grupo != "Control") %>% droplevels()
  
  # Kruskal-Wallis global
  kw_res <- kruskal.test(Value ~ Grupo, data = stats_data)
  
  # Dunn Post-Hoc con ajuste Benjamini-Hochberg (BH)
  capture.output(dunn_res <- dunn.test(stats_data$Value, stats_data$Grupo, method = "bh", kw = FALSE, table = FALSE))
  
  dunn_df <- data.frame(
    Dia = day_name,
    Comparisons = dunn_res$comparisons,
    Z_stat = dunn_res$Z,
    P_raw = dunn_res$P,
    P_adj_BH = dunn_res$P.adjusted
  )
  
  dunn_results_list[[day_name]] <<- dunn_df
  
  max_y <- max(plot_data$Value, na.rm = TRUE)
  min_y <- min(plot_data$Value, na.rm = TRUE)
  range_y <- max_y - min_y
  if (range_y == 0) range_y <- 1
  
  comp_metadata <- list(
    list(g1 = "25ug",  x1 = 2, x2 = 5, row = 1, pair_pattern = "25ug - Comercial|Comercial - 25ug"),
    list(g1 = "50ug",  x1 = 3, x2 = 5, row = 2, pair_pattern = "50ug - Comercial|Comercial - 50ug"),
    list(g1 = "100ug", x1 = 4, x2 = 5, row = 3, pair_pattern = "100ug - Comercial|Comercial - 100ug")
  )
  
  lines_and_labels <- map_df(comp_metadata, function(m) {
    p_val <- dunn_df %>% 
      filter(grepl(m$pair_pattern, Comparisons)) %>% 
      pull(P_adj_BH)
    
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
  
  p <- ggplot(plot_data, aes(x = Grupo, y = Value, color = Grupo)) +
    geom_boxplot(aes(fill = Grupo), width = 0.5, outlier.shape = NA, alpha = 0.14, linewidth = 0.4) +
    geom_jitter(width = 0.1, height = 0, size = 1.3, alpha = 0.8) +
    
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
    coord_cartesian(ylim = c(min_y - (range_y * 0.05), ylim_max)) +
    
    labs(title = meta$title, x = NULL, y = meta$ylabel) +
    theme_nature(11) +
    theme(plot.margin = margin(t = 5, r = 8, b = 5, l = 5))
  
  return(p)
}

plot_list <- imap(metric_metadata, ~ generate_endpoint_titer_boxplot(.x, .y))
fig_titulos <- patchwork::wrap_plots(plot_list, nrow = 1) + 
  patchwork::plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))

save_nature(fig_titulos, "Figure_Titulo_end_point_BH", width = 7, height = 3.8)


# ============================================================
# 6. FIGURA 2: COMPARACIÓN DE TIEMPOS (MANN-WHITNEY U + AJUSTE BH)
# ============================================================

generate_time_comparison_boxplot <- function(datos) {
  grupos_eval <- levels(datos$Grupo)
  
  raw_stats <- map_df(seq_along(grupos_eval), function(i) {
    grp <- grupos_eval[i]
    sub_df <- datos %>% filter(Grupo == grp)
    
    if (sd(sub_df$Log2_Titulo, na.rm = TRUE) == 0 || length(unique(sub_df$Tiempo_Num)) < 2) {
      p_raw <- 1.0
      w_stat <- NA
    } else {
      wilcox_res <- wilcox.test(Log2_Titulo ~ Tiempo_Factor, data = sub_df, exact = FALSE)
      p_raw <- wilcox_res$p.value
      w_stat <- as.numeric(wilcox_res$statistic)
    }
    
    max_grp_y <- max(sub_df$Log2_Titulo, na.rm = TRUE)
    
    data.frame(
      Grupo = grp,
      x1 = i - 0.18,
      x2 = i + 0.18,
      y_bar = max_grp_y + 0.8,
      x_text = i,
      y_text = max_grp_y + 1.2,
      W_statistic = w_stat,
      p_raw = p_raw
    )
  })
  
  stats_list <- raw_stats %>%
    mutate(
      p_adj = p.adjust(p_raw, method = "BH"),
      p_label = case_when(
        p_adj < 0.001 ~ "p < 0.001",
        p_adj >= 0.05 ~ "ns",
        TRUE ~ paste0("p = ", sprintf("%.3f", p_adj))
      )
    )
  
  max_y_total <- max(stats_list$y_text) + 0.5
  min_y_total <- min(datos$Log2_Titulo, na.rm = TRUE)
  
  p <- ggplot(datos, aes(x = Grupo, y = Log2_Titulo, fill = Tiempo_Factor)) +
    geom_boxplot(
      position = position_dodge(width = 0.75), 
      width = 0.6, outlier.shape = NA, linewidth = 0.4, 
      color = "grey20", alpha = 0.85
    ) +
    geom_point(
      aes(group = Tiempo_Factor, color = Tiempo_Factor), 
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), 
      size = 1.2, alpha = 0.9, show.legend = FALSE
    ) +
    geom_segment(data = stats_list, aes(x = x1, xend = x2, y = y_bar, yend = y_bar), color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_segment(data = stats_list, aes(x = x1, xend = x1, y = y_bar, yend = y_bar - 0.2), color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_segment(data = stats_list, aes(x = x2, xend = x2, y = y_bar, yend = y_bar - 0.2), color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_text(data = stats_list, aes(x = x_text, y = y_text, label = p_label), color = "black", size = 2.3, vjust = 0, fontface = "plain", inherit.aes = FALSE) +
    
    scale_fill_manual(name = "Time Point", values = c("Día 28" = "#80B1D3", "Día 42" = "#2B5C8F")) +
    scale_color_manual(values = c("Día 28" = "#2C3E50", "Día 42" = "#1A252F")) +
    coord_cartesian(ylim = c(min_y_total - 0.5, max_y_total)) +
    labs(
      title = "Temporal Comparison of Antibody Titers Across Doses",
      x = NULL,
      y = expression("Endpoint Antibody Titer (" * log[2] * ")")
    ) +
    theme_nature(11) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 7.5, face = "bold"),
      legend.text = element_text(size = 7),
      plot.margin = margin(t = 8, r = 8, b = 5, l = 5)
    )
  
  return(list(plot = p, stats = stats_list))
}

fig_time_res <- generate_time_comparison_boxplot(datos_titulos)
save_nature(fig_time_res$plot, "Figure_Time_Comparison_By_Dose_Boxplot_BH", width = 7, height = 3.8)


# ============================================================
# 7. EXPORTACIÓN DE RESULTADOS A EXCEL
# ============================================================

wb_titulos <- createWorkbook()

# Hoja 1: Resumen de Títulos
resumen_titulos <- datos_titulos %>%
  group_by(Grupo, Tiempo_Factor) %>%
  summarise(
    n = n(),
    Median_Log2 = median(Log2_Titulo, na.rm = TRUE),
    IQR_Log2 = IQR(Log2_Titulo, na.rm = TRUE),
    Mean_Log2 = mean(Log2_Titulo, na.rm = TRUE),
    SD_Log2 = sd(Log2_Titulo, na.rm = TRUE),
    .groups = "drop"
  )
addWorksheet(wb_titulos, "Resumen_Titulos")
writeData(wb_titulos, "Resumen_Titulos", resumen_titulos)

# Hoja 2: Comparaciones Post-Hoc Dunn (Entre Grupos)
df_dunn_all <- bind_rows(dunn_results_list)
addWorksheet(wb_titulos, "Dunn_PostHoc_BH")
writeData(wb_titulos, "Dunn_PostHoc_BH", df_dunn_all)

# Hoja 3: Comparaciones Temporales Mann-Whitney U (Entre Días)
df_wilcox <- fig_time_res$stats %>%
  select(Grupo, W_statistic, p_raw, p_adj_BH = p_adj, p_label)
addWorksheet(wb_titulos, "Mann_Whitney_Temporal_BH")
writeData(wb_titulos, "Mann_Whitney_Temporal_BH", df_wilcox)

saveWorkbook(wb_titulos, file.path(tables_dir, "Resultados_Titulos_No_Parametrico.xlsx"), overwrite = TRUE)

cat("\n[✓] Análisis de Títulos completado exitosamente.")
cat("\n Archivos guardados en:")
cat("\n - Tablas: ", file.path(tables_dir, "Resultados_Titulos_No_Parametrico.xlsx"))
cat("\n - Figuras: ", figures_dir, "\n")
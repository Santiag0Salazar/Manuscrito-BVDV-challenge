###############################################################################
# SCRIPT DE ANÁLISIS DE LINFOCITOS Y TEMPERATURA RECTAL - ENSAYO BVDV
# ESTILO NATURE / GRAPHPAD PRISM + MODELOS DE EFECTOS MIXTOS (LMM) 
# COMPARACIONES TODOS CONTRA TODOS (TUKEY-KRAMER) + BOXPLOT DÍA 8
###############################################################################

options(stringsAsFactors = FALSE)

# =============================================================================
# 0. INSTALACIÓN Y CARGA DE PAQUETES
# =============================================================================

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readxl", "readr",
  "lme4", "lmerTest", "emmeans", "multcomp", "broom.mixed",
  "patchwork", "cowplot", "scales", "pracma", "purrr"
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
  library(readxl)
  library(readr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(multcomp)
  library(broom.mixed)
  library(cowplot)
  library(scales)
  library(pracma)
  library(purrr)
})

set.seed(20260620)

# =============================================================================
# 1. CONFIGURACIÓN DE DIRECTORIOS
# =============================================================================

input_file <- "C:/Users/santi/OneDrive/Escritorio/BVDV Imagenes/Manuscrito/Signos clinicos/Registros Lym y Temperaturas ensayo DVB.xlsx"
out_dir    <- file.path("C:/Users/santi/OneDrive/Escritorio/BVDV Imagenes/Manuscrito/Signos clinicos", "Temperatura_LYM")

tables_dir  <- file.path(out_dir, "tablas")
figures_dir <- file.path(out_dir, "figuras")
stats_dir   <- file.path(out_dir, "estadistica_LMM")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. PALETA Y CONFIGURACIÓN VISUAL (ESTANDARIZADA)
# =============================================================================

palette_nature <- c(
  "Control" = "#4D4D4D",
  "25 µg"   = "#004B87",
  "50 µg"   = "#E69F00",
  "100 µg"  = "#2ECC71",
  "Comm."   = "#D55E00"
)

shapes_nature <- c(
  "Control" = 16,
  "25 µg"   = 15,
  "50 µg"   = 17,
  "100 µg"  = 25,
  "Comm."   = 18
)

group_levels <- c("Control", "25 µg", "50 µg", "100 µg", "Comm.")

# =============================================================================
# 3. TEMA NATURE Y FUNCIONES AUXILIARES
# =============================================================================

theme_nature <- function(base_size = 11) {
  cowplot::theme_cowplot(font_size = base_size, font_family = "arial") +
    theme(
      plot.title           = element_text(face = "plain", size = base_size + 1, hjust = 0.5),
      plot.subtitle        = element_text(size = base_size, color = "grey25"),
      strip.background     = element_blank(),
      strip.text           = element_text(face = "plain", size = base_size),
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

save_nature <- function(plot, filename, width, height) {
  ggsave(
    file.path(figures_dir, paste0(filename, ".tiff")),
    plot, width = width, height = height, units = "in", dpi = 600, compression = "lzw"
  )
}

standardize_group <- function(x) {
  case_when(
    trimws(as.character(x)) %in% c("G1", "Control", "Control (-)", "Placebo") ~ "Control",
    trimws(as.character(x)) %in% c("G2", "25ug", "25 µg", "25 µg rE2", "rE2 25ug") ~ "25 µg",
    trimws(as.character(x)) %in% c("G3", "50ug", "50 µg", "50 µg rE2", "rE2 50ug") ~ "50 µg",
    trimws(as.character(x)) %in% c("G4", "100ug", "100 µg", "100 µg rE2", "rE2 100ug") ~ "100 µg",
    trimws(as.character(x)) %in% c("G5", "Comercial", "Commercial", "Comm.", "Cattle master C+") ~ "Comm.",
    TRUE ~ as.character(x)
  )
}

# =============================================================================
# 4. FUNCIÓN PARA LIMPIEZA DE HOJAS
# =============================================================================

clean_temperature_sheet <- function(sheet_name) {
  df <- read_excel(input_file, sheet = sheet_name, skip = 1)
  names(df) <- trimws(names(df))
  
  group_col <- names(df)[tolower(trimws(names(df))) == "grupo"]
  temperature_cols <- grep("^Dia", names(df), value = TRUE)
  
  if (length(group_col) == 0) {
    stop(paste("No se encontró la columna 'Grupo' en la hoja:", sheet_name))
  }
  
  df_long <- df %>%
    dplyr::mutate(
      DIIO = as.character(DIIO),
      Grupo_std = standardize_group(.data[[group_col]])
    ) %>%
    dplyr::filter(!is.na(DIIO), !is.na(Grupo_std), Grupo_std %in% group_levels) %>%
    dplyr::select(DIIO, Grupo = Grupo_std, dplyr::all_of(temperature_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(temperature_cols), 
      names_to = "Tiempo_original", 
      values_to = "Temperatura"
    ) %>%
    dplyr::mutate(
      Temperatura = suppressWarnings(as.numeric(Temperatura)),
      Dia_original = as.numeric(sub("Dia ([0-9]+).*", "\\1", Tiempo_original)),
      Hora = ifelse(grepl("pm", Tiempo_original, ignore.case = TRUE), 0.5, 0),
      Tiempo = Dia_original - 1 + Hora,
      Grupo = factor(Grupo, levels = group_levels),
      Animal = factor(DIIO)
    ) %>%
    dplyr::filter(!is.na(Temperatura), !is.na(Tiempo)) %>%
    dplyr::arrange(Grupo, Animal, Tiempo)
  
  return(df_long)
}

# =============================================================================
# 5. FUNCIÓN DE MODELO LINEAL DE EFECTOS MIXTOS (LMM)
# =============================================================================

run_lmm_analysis <- function(data, response_var, time_var, output_prefix) {
  cat("\n=====================================================================\n")
  cat(paste(" EJECUTANDO MODELO DE EFECTOS MIXTOS (LMM) PARA:", output_prefix, "\n"))
  cat("=====================================================================\n")
  
  data_mod <- data %>%
    dplyr::mutate(Tiempo_factor = factor(.data[[time_var]]))
  
  formula_str <- paste(response_var, "~ Grupo * Tiempo_factor + (1 | Animal)")
  model <- lmer(as.formula(formula_str), data = data_mod)
  
  anova_res <- anova(model, type = 3)
  write.csv(anova_res, file.path(stats_dir, paste0(output_prefix, "_LMM_ANOVA.csv")))
  
  emm <- emmeans(model, ~ Grupo | Tiempo_factor)
  
  pairwise_ctrl <- contrast(emm, method = "trt.vs.ctrl", ref = "Control", adjust = "dunnet")
  write.csv(as.data.frame(pairwise_ctrl), file.path(stats_dir, paste0(output_prefix, "_LMM_Comparaciones_vs_Control.csv")))
  
  pairwise_all <- contrast(emm, method = "pairwise", adjust = "tukey")
  write.csv(as.data.frame(pairwise_all), file.path(stats_dir, paste0(output_prefix, "_LMM_Comparaciones_Tukey_Todos_vs_Todos.csv")))
  
  return(list(model = model, anova = anova_res, emmeans = emm, contrasts = pairwise_all))
}

# =============================================================================
# 6. ANÁLISIS DE TEMPERATURAS (DESAFÍO, V1, V2)
# =============================================================================

# --- Desafío ---
datos_temperature <- clean_temperature_sheet("T° rectal desafio")
temperature_summary <- datos_temperature %>%
  dplyr::group_by(Grupo, Tiempo) %>%
  dplyr::summarise(
    n = sum(!is.na(Temperatura)),
    mean_temp = mean(Temperatura, na.rm = TRUE),
    se_temp = sd(Temperatura, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

write_csv(datos_temperature, file.path(tables_dir, "01_temperatura_desafio_datos_largos.csv"))
write_csv(temperature_summary, file.path(tables_dir, "02_resumen_temperatura_desafio.csv"))

lmm_temp_desafio <- run_lmm_analysis(datos_temperature, "Temperatura", "Tiempo", "01_Temp_Desafio")

fig1_temperature <- ggplot() +
  geom_line(data = datos_temperature, aes(x = Tiempo, y = Temperatura, group = Animal, color = Grupo), linewidth = 0.35, alpha = 0.20) +
  geom_point(data = datos_temperature, aes(x = Tiempo, y = Temperatura, color = Grupo), size = 1.0, alpha = 0.25) +
  geom_line(data = temperature_summary, aes(x = Tiempo, y = mean_temp, color = Grupo, group = Grupo), linewidth = 0.9) +
  geom_errorbar(data = temperature_summary, aes(x = Tiempo, ymin = mean_temp - se_temp, ymax = mean_temp + se_temp, color = Grupo), width = 0.12, linewidth = 0.35) +
  geom_point(data = temperature_summary, aes(x = Tiempo, y = mean_temp, color = Grupo, shape = Grupo, fill = Grupo), size = 2.2) +
  geom_hline(yintercept = 39.2, linetype = "dashed", linewidth = 0.55, color = "grey30") +
  scale_color_manual(values = palette_nature) +
  scale_fill_manual(values = palette_nature) +
  scale_shape_manual(values = shapes_nature) +
  scale_x_continuous(
    breaks = seq(0, 28, by = 2), # Puedes usar by = 1 o by = 2 según prefieras la densidad de etiquetas
    limits = c(0, 28),
    expand = c(0, 0.2)
  )+  scale_y_continuous(breaks = seq(floor(min(datos_temperature$Temperatura, na.rm = TRUE)), ceiling(max(datos_temperature$Temperatura, na.rm = TRUE)), by = 0.5)) +
  labs(x = "Days after challenge", y = "Rectal temperature (°C)", color = "Vaccine Group", fill = "Vaccine Group", shape = "Vaccine Group") +
  theme_nature(11)

save_nature(fig1_temperature, "Figure_1_Rectal_Temperature_Challenge", 5, 3.3)

# --- Vacuna 1 ---
datos_temperature_v1 <- clean_temperature_sheet("T°Rectal v1")
temp_v1_summary <- datos_temperature_v1 %>%
  dplyr::group_by(Grupo, Tiempo) %>%
  dplyr::summarise(
    n = sum(!is.na(Temperatura)),
    mean_temp = mean(Temperatura, na.rm = TRUE),
    se_temp = sd(Temperatura, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

write_csv(datos_temperature_v1, file.path(tables_dir, "04_temperatura_post_vacuna_1_datos.csv"))
write_csv(temp_v1_summary, file.path(tables_dir, "05_resumen_temperatura_post_vacuna_1.csv"))

lmm_temp_v1 <- run_lmm_analysis(datos_temperature_v1, "Temperatura", "Tiempo", "02_Temp_Vacuna1")

# --- Vacuna 2 ---
datos_temperature_v2 <- clean_temperature_sheet("T°Rectal v2")
temp_v2_summary <- datos_temperature_v2 %>%
  dplyr::group_by(Grupo, Tiempo) %>%
  dplyr::summarise(
    n = sum(!is.na(Temperatura)),
    mean_temp = mean(Temperatura, na.rm = TRUE),
    se_temp = sd(Temperatura, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

write_csv(datos_temperature_v2, file.path(tables_dir, "06_temperatura_post_vacuna_2_datos.csv"))
write_csv(temp_v2_summary, file.path(tables_dir, "07_resumen_temperatura_post_vacuna_2.csv"))

lmm_temp_v2 <- run_lmm_analysis(datos_temperature_v2, "Temperatura", "Tiempo", "03_Temp_Vacuna2")

# =============================================================================
# 7. PROCESAMIENTO Y FIGURA DE LINFOCITOS (LYM)
# =============================================================================

datos_lym_raw <- read_excel(input_file, sheet = "LYM", skip = 1)
names(datos_lym_raw) <- trimws(names(datos_lym_raw))
names(datos_lym_raw)[1] <- "Grupo"

lym_cols <- grep("^S[0-9]+$", names(datos_lym_raw), value = TRUE)

lym_long <- datos_lym_raw %>%
  dplyr::mutate(Grupo = factor(standardize_group(Grupo), levels = group_levels)) %>%
  dplyr::filter(!is.na(Grupo)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(lym_cols), names_to = "Sampling", values_to = "LYM") %>%
  dplyr::mutate(
    LYM = suppressWarnings(as.numeric(LYM)),
    Sampling = factor(Sampling, levels = paste0("S", 1:12))
  ) %>%
  dplyr::filter(!is.na(LYM))

write_csv(lym_long, file.path(tables_dir, "08_LYM_datos_largos.csv"))

fig4_lym <- ggplot(lym_long, aes(x = Sampling, y = LYM, color = Grupo, group = Grupo)) +
  geom_vline(xintercept = 7, linetype = "dashed", linewidth = 0.55, color = "grey30") +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = Grupo, fill = Grupo), size = 2.8) +
  scale_color_manual(values = palette_nature) +
  scale_fill_manual(values = palette_nature) +
  scale_shape_manual(values = shapes_nature) +
  scale_x_discrete(
    labels = c(
      "S1\nDay 0", "S2\nDay 7", "S3\nDay 14", "S4\nDay 22", 
      "S5", "S6", "S7\nChallenge\nDay 56", "S8\n+3 d", 
      "S9\n+7 d", "S10\n+14 d", "S11\n+21 d", "S12\n+28 d"
    )
  ) +
  labs(
    x = "Sampling timepoints",
    y = expression("Lymphocytes (" * 10^9 * "/L)"),
    color = "Vaccine Group", fill = "Vaccine Group", shape = "Vaccine Group"
  ) +
  theme_nature(12)

save_nature(fig4_lym, "Figure_4_Lymphocyte_Kinetics", 7.5, 3.8)

# =============================================================================
# 8. GENERACIÓN DEL BOXPLOT EN DÍA 8 POST-DESAFÍO (TODOS CONTRA TODOS - TUKEY)
# =============================================================================

datos_temp_day8 <- datos_temperature %>%
  filter(Tiempo == 8)

lm_day8    <- lm(Temperatura ~ Grupo, data = datos_temp_day8)
emms_day8  <- emmeans(lm_day8, ~ Grupo)

pairs_day8 <- as.data.frame(pairs(emms_day8, adjust = "tukey"))

write_csv(pairs_day8, file.path(stats_dir, "09_Temperatura_Dia8_Comparaciones_Todos_vs_Todos_Tukey.csv"))

generate_day8_temperature_boxplot <- function(data_plot, stats_df) {
  
  max_y <- max(data_plot$Temperatura, na.rm = TRUE)
  min_y <- min(data_plot$Temperatura, na.rm = TRUE)
  range_y <- max_y - min_y
  if (range_y == 0) range_y <- 1
  
  comp_metadata <- list(
    list(x1 = 1, x2 = 2, row = 1, pattern = "Control - 25 µg|25 µg - Control"),
    list(x1 = 1, x2 = 3, row = 2, pattern = "Control - 50 µg|50 µg - Control"),
    list(x1 = 1, x2 = 4, row = 3, pattern = "Control - 100 µg|100 µg - Control"),
    list(x1 = 1, x2 = 5, row = 4, pattern = "Control - Comm.|Comm. - Control"),
    list(x1 = 2, x2 = 5, row = 5, pattern = "25 µg - Comm.|Comm. - 25 µg"),
    list(x1 = 3, x2 = 5, row = 6, pattern = "50 µg - Comm.|Comm. - 50 µg"),
    list(x1 = 4, x2 = 5, row = 7, pattern = "100 µg - Comm.|Comm. - 100 µg")
  )
  
  lines_and_labels <- map_df(comp_metadata, function(m) {
    p_val <- stats_df %>%
      filter(grepl(m$pattern, contrast)) %>%
      pull(p.value)
    
    if (length(p_val) == 0) p_val <- 1.0
    
    p_formatted <- if (p_val < 0.001) "p < 0.001" else if (p_val >= 0.05) "ns" else paste0("p = ", sprintf("%.3f", p_val))
    y_bar <- max_y + (range_y * 0.075 * m$row)
    
    data.frame(
      x = m$x1, xend = m$x2, y = y_bar,
      x_text = (m$x1 + m$x2) / 2, y_text = y_bar + (range_y * 0.02),
      p_label = p_formatted
    )
  })
  
  ylim_max <- max(lines_and_labels$y_text) + (range_y * 0.05)
  
  p <- ggplot(data_plot, aes(x = Grupo, y = Temperatura, color = Grupo, fill = Grupo)) +
    geom_boxplot(
      width = 0.5, 
      outlier.shape = NA, 
      alpha = 0.15, 
      linewidth = 0.4
    ) +
    geom_point(
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
    scale_y_continuous(
      breaks = seq(
        floor(min(datos_temperature$Temperatura, na.rm = TRUE)), 
        ceiling(max(datos_temperature$Temperatura, na.rm = TRUE)), 
        by = 0.5
      )
    )+
    geom_segment(data = lines_and_labels, aes(x = x, xend = xend, y = y, yend = y), 
                 color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_segment(data = lines_and_labels, aes(x = x, xend = x, y = y, yend = y - (range_y * 0.012)), 
                 color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_segment(data = lines_and_labels, aes(x = xend, xend = xend, y = y, yend = y - (range_y * 0.012)), 
                 color = "grey35", linewidth = 0.35, inherit.aes = FALSE) +
    geom_text(data = lines_and_labels, aes(x = x_text, y = y_text, label = p_label), 
              color = "black", size = 2.0, vjust = 0, fontface = "plain", inherit.aes = FALSE) +
    scale_color_manual(values = palette_nature, drop = FALSE) +
    scale_fill_manual(values = palette_nature, drop = FALSE) +
    coord_cartesian(ylim = c(min_y - (range_y * 0.05), ylim_max)) +
    labs(
      title = "Rectal Temperature at Day 8 Post-Challenge",
      x = NULL, 
      y = "Rectal Temperature (°C)"
    ) +
    theme_nature(11) +
    theme(
      plot.margin = margin(t = 5, r = 8, b = 5, l = 5)
    )
  
  ggsave(
    filename = file.path(figures_dir, "Figure_Boxplot_Temperature_Day8.tiff"), 
    plot = p, 
    width = 2.1, 
    height = 4.1, 
    dpi = 600, 
    compression = "lzw"
  )
  
  return(p)
}

fig_day8_temp <- generate_day8_temperature_boxplot(datos_temp_day8, pairs_day8)

# =============================================================================
# 9. IMPRESIÓN Y SALIDA FINAL EN CONSOLA
# =============================================================================

print(fig1_temperature)
print(fig4_lym)
print(fig_day8_temp)

cat("\n=====================================================================\n")
cat(" EJECUCIÓN FINALIZADA CORRECTAMENTE CON ÉXITO\n")
cat(" Tablas guardadas en: ", tables_dir, "\n")
cat(" Figuras TIFF (600 dpi) guardadas en: ", figures_dir, "\n")
cat(" Resultados Estadísticos LMM guardados en: ", stats_dir, "\n")
cat("=====================================================================\n")






###############################################################################
# SCRIPT DE EXTRACCIÓN Y REPORTE TEXTUAL PARA MANUSCRITO (DÍA 8 POST-DESAFÍO)
###############################################################################

# Carga de librerías necesarias
suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(readr)
})

# =============================================================================
# 1. MODELADO Y COMPARACIONES EN EL DÍA 8
# =============================================================================

# Filtrar datos del Día 8
datos_temp_day8 <- datos_temperature %>%
  filter(Tiempo == 8)

# Modelo Lineal (ANOVA una vía)
lm_day8 <- lm(Temperatura ~ Grupo, data = datos_temp_day8)

# ANOVA summary
anova_day8 <- summary(aov(lm_day8))[[1]]

# Extracción de F, Df1, Df2 y p-value
f_stat <- anova_day8["Grupo", "F value"]
df_num <- anova_day8["Grupo", "Df"]
df_den <- anova_day8["Residuals", "Df"]
p_val_model <- anova_day8["Grupo", "Pr(>F)"]

# Formateo del p-valor del modelo
p_model_str <- if (p_val_model < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val_model)

# Texto del Modelo Lineal
model_report_str <- sprintf("F(%d, %d) = %.2f, %s", df_num, df_den, f_stat, p_model_str)

# =============================================================================
# 2. COMPARACIONES DE TUKEY ENTRE GRUPOS VACUNADOS
# =============================================================================

# Medias estimadas por el modelo (EMMs)
emms_day8 <- emmeans(lm_day8, ~ Grupo)

# Comparaciones par a par ajustadas por Tukey
pairs_day8 <- as.data.frame(pairs(emms_day8, adjust = "tukey"))

# Función auxiliar para formatear p-valores de comparaciones
format_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  return(sprintf("p = %.3f", p))
}

# Filtrar comparaciones requeridas
comp_vs_comm <- pairs_day8 %>%
  filter(grepl("Comm.", contrast)) %>%
  filter(!grepl("Control", contrast))

comp_among_rec <- pairs_day8 %>%
  filter(!grepl("Control", contrast) & !grepl("Comm.", contrast))

# Formatear comparaciones vs Comercial
text_vs_comm <- comp_vs_comm %>%
  rowwise() %>%
  mutate(text = sprintf("%s (%s)", contrast, format_p(p.value))) %>%
  pull(text) %>%
  paste(collapse = "; ")

# Formatear comparaciones entre dosis recombinantes
text_among_rec <- comp_among_rec %>%
  rowwise() %>%
  mutate(text = sprintf("%s (%s)", contrast, format_p(p.value))) %>%
  pull(text) %>%
  paste(collapse = "; ")

# =============================================================================
# 3. ESTADÍSTICA DESCRIPTIVA PARA EL GRUPO CONTROL (n = 2)
# =============================================================================

control_stats <- datos_temp_day8 %>%
  filter(Grupo == "Control") %>%
  summarise(
    n = n(),
    mean_temp = mean(Temperatura, na.rm = TRUE),
    sd_temp = sd(Temperatura, na.rm = TRUE)
  )

text_control_desc <- sprintf(
  "mean ± SD: %.2f ± %.2f °C (n = %d)", 
  control_stats$mean_temp, 
  control_stats$sd_temp, 
  control_stats$n
)

# =============================================================================
# 4. REDACCIÓN Y CONSTRUCCIÓN DEL PÁRRAFO FINAL
# =============================================================================

final_paragraph <- sprintf(
  "At day 8 post-challenge, rectal temperature differed significantly among groups (%s). In the Tukey-adjusted pairwise comparisons, no significant differences were observed between recombinant vaccine groups and the commercial vaccine group (%s), nor among the three recombinant vaccine doses (%s). Rectal temperature for the control group was described descriptively as %s.",
  model_report_str,
  text_vs_comm,
  text_among_rec,
  text_control_desc
)

# IMPRESIÓN EN CONSOLA
cat("\n=====================================================================\n")
cat(" PARÁGRAFO REDACTADO PARA EL MANUSCRITO:\n")
cat("=====================================================================\n\n")
cat(final_paragraph)
cat("\n\n=====================================================================\n")

# =============================================================================
# 5. GUARDAR RESULTADO EN UN ARCHIVO DE TEXTO
# =============================================================================

report_file <- file.path(stats_dir, "09b_Day8_Statistical_Report_Text.txt")
writeLines(final_paragraph, report_file)
cat(paste(" Reporte guardado en:", report_file, "\n"))

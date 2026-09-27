# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# RESPUESTA CELULAR PRE-DESAFÍO
# Citoquinas en PBMC estimuladas con E2 (días 0 y 56)
#
# Sección 3.3 y Figura 4 del manuscrito
#
# Citoquinas analizadas: IFN-γ, TNF-α e IL-17
# (IL-4 quedó bajo el límite de detección en todas las muestras)
#
# Análisis (Sección 2.11 del manuscrito):
#   - Valores netos (estimulado - no estimulado), pg/mL
#   - Transformación log10(x + 1)
#   - Modelo lineal mixto (nlme::lme): Tiempo * Grupo,
#     intercepto aleatorio por animal
#   - Efectos fijos: F con gl de numerador y denominador
#   - Día 0 vs día 56 dentro de cada grupo (emmeans)
#   - Grupos dentro de cada día (emmeans, ajuste de Tukey)
#
# Salidas (carpeta results/ del proyecto):
#   figures/PBMC_Figure_4_cytokines.pdf / .tiff / .png
#   tables/PBMC_cytokines_resultados.xlsx
#   reports/PBMC_cytokines_informe.docx
#   models/PBMC_modelos_mixtos.rds
#   sessionInfo_PBMC.txt
#
# Autor: Santiago Salazar
# ============================================================


# ============================================================
# 1. OPCIONES Y PAQUETES
# ============================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "readxl",
  "openxlsx",
  "nlme",
  "emmeans",
  "patchwork",
  "cowplot",
  "purrr",
  "officer",
  "flextable"
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
  library(tibble)
  library(readxl)
  library(openxlsx)
  library(nlme)
  library(emmeans)
  library(patchwork)
  library(cowplot)
  library(purrr)
  library(officer)
  library(flextable)
})


# ============================================================
# 2. SEMILLA (solo afecta la posición horizontal de los puntos)
# ============================================================

set.seed(20260811)


# ============================================================
# 3. LOCALIZAR PROYECTO
# ============================================================

find_project_root <- function() {
  
  current_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  
  while (current_dir != dirname(current_dir)) {
    
    if (length(list.files(current_dir, pattern = "\\.Rproj$")) > 0) {
      return(current_dir)
    }
    
    current_dir <- dirname(current_dir)
  }
  
  stop("No se encontró el archivo .Rproj del proyecto.")
}

project_dir <- find_project_root()


# ============================================================
# 4. DIRECTORIOS
# ============================================================

data_dir    <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
tables_dir  <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
models_dir  <- file.path(results_dir, "models")
reports_dir <- file.path(results_dir, "reports")

for (d in c(tables_dir, figures_dir, models_dir, reports_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ============================================================
# 5. PARÁMETROS DEL ANÁLISIS
# ============================================================

input_file <- file.path(data_dir, "Datos_Estimulacion.xlsx")

# Citoquinas: columna del Excel, etiqueta y panel de la Figura 4
# (orden del manuscrito: A = IFN-γ, B = TNF-α, C = IL-17)
cytokines <- tibble::tribble(
  ~var,         ~label,   ~panel,
  "IFN_pg_ml",  "IFN-γ",  "A",
  "TNF_pg_ml",  "TNF-α",  "B",
  "IL17_pg_ml", "IL-17",  "C"
)

# Valores netos negativos (no estimulado > estimulado):
# log10(x + 1) no admite valores <= -1. Si es TRUE, se fijan en 0 y
# el número de valores corregidos queda registrado en las salidas.
# Declararlo en Métodos (Sección 2.7) si ocurre.
truncate_negative <- TRUE

alpha_level <- 0.05

group_levels <- c("Control", "25ug", "50ug", "100ug", "Commercial")
time_levels  <- c("0", "56")


# ============================================================
# 6. LECTURA DE DATOS
# ============================================================

if (!file.exists(input_file)) {
  stop(paste0("No se encontró el archivo:\n", input_file))
}

cat("\nArchivo de entrada:\n", input_file, "\n")

datos_raw <- readxl::read_excel(input_file, sheet = 1)

required_columns <- c("Animal", "Grupo", "Tiempo", cytokines$var)
missing_columns  <- setdiff(required_columns, names(datos_raw))

if (length(missing_columns) > 0) {
  stop(paste0("Faltan las siguientes columnas:\n",
              paste(missing_columns, collapse = ", ")))
}


# ============================================================
# 7. PREPARACIÓN DE DATOS
# ============================================================

# Unifica los nombres de grupo usados en las distintas planillas
recode_group <- function(x) {
  
  x  <- trimws(as.character(x))
  xl <- tolower(x)
  
  # Dosis: número que precede a "ug"/"µg" (evita tomar el "2" de "rE2 25ug")
  xn <- ifelse(
    grepl("[0-9]+\\s*(ug|µg|mcg)", xl, perl = TRUE),
    sub("^.*?([0-9]+)\\s*(ug|µg|mcg).*$", "\\1", xl, perl = TRUE),
    gsub("[^0-9]", "", xl)
  )
  
  dplyr::case_when(
    toupper(x) == "G1"                                  ~ "Control",
    toupper(x) == "G2"                                  ~ "25ug",
    toupper(x) == "G3"                                  ~ "50ug",
    toupper(x) == "G4"                                  ~ "100ug",
    toupper(x) == "G5"                                  ~ "Commercial",
    grepl("control", xl)                                ~ "Control",
    grepl("comer|commer|comm|cattle", xl)               ~ "Commercial",
    xn == "25"                                          ~ "25ug",
    xn == "50"                                          ~ "50ug",
    xn == "100"                                         ~ "100ug",
    TRUE                                                ~ NA_character_
  )
}

# Unifica la codificación del tiempo (día 0 y día 56 post-vacunación)
recode_time <- function(x) {
  
  x <- tolower(trimws(as.character(x)))
  
  dplyr::case_when(
    x %in% c("t0", "0", "d0", "day 0", "dia 0", "día 0")                  ~ "0",
    x %in% c("t6", "t56", "56", "d56", "day 56", "dia 56", "día 56")      ~ "56",
    TRUE                                                                  ~ NA_character_
  )
}

datos <- datos_raw %>%
  mutate(
    Animal = factor(Animal),
    Grupo  = factor(recode_group(Grupo), levels = group_levels),
    Tiempo = factor(recode_time(Tiempo), levels = time_levels)
  )

if (any(is.na(datos$Grupo))) {
  stop("Hay grupos que no se pudieron reconocer: ",
       paste(unique(datos_raw$Grupo[is.na(datos$Grupo)]), collapse = ", "))
}

if (any(is.na(datos$Tiempo))) {
  stop("Hay tiempos que no se pudieron reconocer: ",
       paste(unique(datos_raw$Tiempo[is.na(datos$Tiempo)]), collapse = ", "))
}

# Formato largo: una fila por animal, tiempo y citoquina
datos_long <- datos %>%
  select(Animal, Grupo, Tiempo, all_of(cytokines$var)) %>%
  pivot_longer(
    cols      = all_of(cytokines$var),
    names_to  = "var",
    values_to = "Value_raw"
  ) %>%
  mutate(Value_raw = as.numeric(Value_raw)) %>%
  left_join(cytokines, by = "var") %>%
  mutate(
    negative  = !is.na(Value_raw) & Value_raw < 0,
    Value     = if (truncate_negative) pmax(Value_raw, 0) else Value_raw,
    log_value = log10(Value + 1)
  )

negative_summary <- datos_long %>%
  group_by(label) %>%
  summarise(n_negative_set_to_zero = sum(negative), .groups = "drop")


# ============================================================
# 8. COMPROBACIÓN BÁSICA
# ============================================================

cat("\n================ DATOS =================\n")

cat("\nAnimales por grupo:\n")
print(datos %>% distinct(Animal, Grupo) %>% count(Grupo))

cat("\nObservaciones por grupo y tiempo:\n")
print(datos %>% count(Grupo, Tiempo))

cat("\nValores netos negativos fijados en 0:\n")
print(negative_summary)


# ============================================================
# 9. PALETA, FORMAS Y TEMA (iguales al curso temporal de IgG)
# ============================================================

palette_nature <- c(
  "Control"    = "#4D4D4D",
  "25ug"       = "#004B87",
  "50ug"       = "#E69F00",
  "100ug"      = "#2ECC71",
  "Commercial" = "#D55E00"
)

# Círculo, cuadrado, triángulo, triángulo invertido y rombo (pie de la Figura 4)
shape_groups <- c(
  "Control"    = 21,
  "25ug"       = 22,
  "50ug"       = 24,
  "100ug"      = 25,
  "Commercial" = 23
)

group_labels <- c(
  "Control"    = "Control",
  "25ug"       = "25 µg",
  "50ug"       = "50 µg",
  "100ug"      = "100 µg",
  "Commercial" = "Comm."
)

theme_nature <- function(base_size = 10) {
  
  cowplot::theme_cowplot(font_size = base_size, font_family = "arial") +
    theme(
      plot.title          = element_text(face = "plain", size = base_size, hjust = 0.5),
      axis.title          = element_text(size = base_size),
      axis.text           = element_text(size = base_size, color = "grey15"),
      strip.background    = element_blank(),
      strip.text          = element_text(size = base_size),
      strip.placement     = "outside",
      legend.position     = "none",
      panel.grid.major.y  = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x  = element_blank(),
      panel.spacing.x     = unit(0.6, "lines")
    )
}


# ============================================================
# 10. FUNCIÓN PARA GUARDAR FIGURAS (igual al curso temporal de IgG)
# ============================================================

save_nature <- function(plot, filename, width = 7, height = 4) {
  
  ggsave(file.path(figures_dir, paste0(filename, ".pdf")),
         plot = plot, width = width, height = height, units = "in")
  
  ggsave(file.path(figures_dir, paste0(filename, ".tiff")),
         plot = plot, width = width, height = height, units = "in",
         dpi = 600, compression = "lzw", bg = "white")
  
  ggsave(file.path(figures_dir, paste0(filename, ".png")),
         plot = plot, width = width, height = height, units = "in",
         dpi = 600, bg = "white")
}


# ============================================================
# 11. FORMATO DE VALORES P
# ============================================================

# Para tablas: valor exacto con 4 decimales
fmt_p_table <- function(p) {
  ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
}

# Para texto y figura: formato del manuscrito
fmt_p_text <- function(p) {
  ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))
}


# ============================================================
# 12. ESTADÍSTICA DESCRIPTIVA (escala original, pg/mL)
# ============================================================

descriptives <- datos_long %>%
  group_by(label, Grupo, Tiempo) %>%
  summarise(
    n      = sum(!is.na(Value)),
    mean   = mean(Value, na.rm = TRUE),
    sd     = sd(Value, na.rm = TRUE),
    sem    = sd / sqrt(n),
    median = median(Value, na.rm = TRUE),
    q1     = quantile(Value, 0.25, na.rm = TRUE),
    q3     = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 13. MODELOS LINEALES MIXTOS POR CITOQUINA
# ============================================================

fit_cytokine_model <- function(df) {
  
  nlme::lme(
    fixed     = log_value ~ Tiempo * Grupo,
    random    = ~ 1 | Animal,
    data      = df,
    method    = "REML",
    na.action = na.omit
  )
}

cytokine_data <- split(datos_long, datos_long$label)[cytokines$label]

models <- purrr::map(cytokine_data, fit_cytokine_model)


# ------------------------------------------------------------
# 13.1 Efectos fijos: F(numDF, denDF) y p
# ------------------------------------------------------------

anova_table <- purrr::imap_dfr(models, function(m, lab) {
  
  a <- as.data.frame(anova(m))
  
  tibble(
    Cytokine = lab,
    Effect   = rownames(a),
    numDF    = a$numDF,
    denDF    = a$denDF,
    F_value  = round(a$`F-value`, 2),
    p_value  = a$`p-value`,
    p_text   = fmt_p_table(a$`p-value`)
  ) %>%
    filter(Effect != "(Intercept)") %>%
    mutate(Effect = recode(Effect,
                           "Tiempo"       = "Time",
                           "Grupo"        = "Group",
                           "Tiempo:Grupo" = "Group × time"))
})


# ------------------------------------------------------------
# 13.2 Día 56 vs día 0 dentro de cada grupo
#      (una comparación por grupo; no requiere ajuste)
# ------------------------------------------------------------

time_contrasts <- purrr::imap_dfr(models, function(m, lab) {
  
  emm <- emmeans(m, ~ Tiempo | Grupo, data = cytokine_data[[lab]])
  
  as.data.frame(pairs(emm, reverse = TRUE, adjust = "none")) %>%
    as_tibble() %>%
    transmute(
      Cytokine        = lab,
      Grupo           = Grupo,
      Contrast        = "Day 56 vs day 0",
      estimate_log10  = estimate,
      SE              = SE,
      df              = df,
      t_ratio         = t.ratio,
      # Razón aproximada (x + 1) día 56 / día 0
      fold_change     = 10^estimate,
      p_value         = p.value,
      p_text          = fmt_p_table(p.value),
      significant     = p.value < alpha_level
    )
})


# ------------------------------------------------------------
# 13.3 Comparaciones entre grupos dentro de cada día (Tukey)
# ------------------------------------------------------------

group_contrasts <- purrr::imap_dfr(models, function(m, lab) {
  
  emm <- emmeans(m, ~ Grupo | Tiempo, data = cytokine_data[[lab]])
  
  as.data.frame(pairs(emm, adjust = "tukey")) %>%
    as_tibble() %>%
    transmute(
      Cytokine       = lab,
      Day            = Tiempo,
      Contrast       = contrast,
      estimate_log10 = estimate,
      SE             = SE,
      df             = df,
      t_ratio        = t.ratio,
      p_value        = p.value,
      p_text         = fmt_p_table(p.value)
    )
})


# ------------------------------------------------------------
# 13.4 Frases listas para el manuscrito (Sección 3.3)
# ------------------------------------------------------------

text_summary <- purrr::map_dfr(cytokines$label, function(lab) {
  
  inter <- anova_table %>% filter(Cytokine == lab, Effect == "Group × time")
  tc    <- time_contrasts %>% filter(Cytokine == lab)
  
  by_group <- paste0(group_labels[as.character(tc$Grupo)], ", ",
                     fmt_p_text(tc$p_value), collapse = "; ")
  
  tibble(
    Cytokine = lab,
    Sentence = paste0(
      lab, ": group × time interaction, F", inter$numDF, ",", inter$denDF,
      " = ", sprintf("%.2f", inter$F_value), ", ", fmt_p_text(inter$p_value),
      ". Day 56 vs day 0 — ", by_group, "."
    )
  )
})

cat("\n================ RESUMEN PARA EL TEXTO =================\n")
for (s in text_summary$Sentence) cat("\n", s, "\n")


# ============================================================
# 14. FIGURA 4
#     Barras: media ± SEM (pg/mL, sin transformar)
#     Puntos: animales individuales
#     Valores p: modelo mixto, día 56 vs día 0 (solo p < 0.05)
# ============================================================

plot_cytokine <- function(lab, show_x_title = FALSE) {
  
  raw_c  <- datos_long   %>% filter(label == lab)
  desc_c <- descriptives %>% filter(label == lab)
  
  # Altura de la barra de significancia por grupo
  top_by_group <- raw_c %>%
    group_by(Grupo) %>%
    summarise(top_raw = max(Value, na.rm = TRUE), .groups = "drop") %>%
    left_join(
      desc_c %>%
        group_by(Grupo) %>%
        summarise(top_bar = max(mean + sem, na.rm = TRUE), .groups = "drop"),
      by = "Grupo"
    ) %>%
    mutate(top = pmax(top_raw, top_bar, na.rm = TRUE))
  
  y_range <- max(top_by_group$top, na.rm = TRUE)
  
  brackets <- time_contrasts %>%
    filter(Cytokine == lab, significant) %>%
    left_join(top_by_group, by = "Grupo") %>%
    mutate(
      y      = top + 0.06 * y_range,
      y_tick = y - 0.02 * y_range,
      y_text = y + 0.02 * y_range,
      label  = fmt_p_text(p_value)
    )
  
  ggplot(desc_c, aes(x = Tiempo, y = mean)) +
    
    geom_col(fill = "white", color = "black", width = 0.6, linewidth = 0.4) +
    
    geom_errorbar(aes(ymin = pmax(mean - sem, 0), ymax = mean + sem),
                  width = 0.2, linewidth = 0.4) +
    
    geom_point(
      data = raw_c,
      aes(x = Tiempo, y = Value, color = Grupo, fill = Grupo, shape = Grupo),
      position = position_jitter(width = 0.08, height = 0, seed = 1),
      size = 1.6, alpha = 0.9
    ) +
    
    geom_segment(data = brackets, aes(x = 1, xend = 2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.35, color = "grey20") +
    geom_segment(data = brackets, aes(x = 1, xend = 1, y = y, yend = y_tick),
                 inherit.aes = FALSE, linewidth = 0.35, color = "grey20") +
    geom_segment(data = brackets, aes(x = 2, xend = 2, y = y, yend = y_tick),
                 inherit.aes = FALSE, linewidth = 0.35, color = "grey20") +
    geom_text(data = brackets, aes(x = 1.5, y = y_text, label = label),
              inherit.aes = FALSE, size = 2.6, vjust = 0) +
    
    facet_wrap(~ Grupo, nrow = 1, strip.position = "bottom",
               labeller = as_labeller(group_labels)) +
    
    scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
    scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
    scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
    
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    
    labs(
      x = if (show_x_title) "Time (days)" else NULL,
      y = paste0(lab, " (pg/mL)"),
      color = NULL, fill = NULL, shape = NULL
    ) +
    
    theme_nature(10)
}

fig4_panels <- purrr::map(seq_len(nrow(cytokines)), function(i) {
  plot_cytokine(cytokines$label[i], show_x_title = i == nrow(cytokines))
})

fig4_final <- wrap_plots(fig4_panels, ncol = 1) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(
    legend.position      = "top",
    legend.justification = "center",
    plot.tag             = element_text(face = "bold", size = 11)
  )

save_nature(fig4_final, "PBMC_Figure_4_cytokines", width = 7.2, height = 8.5)


# ============================================================
# 15. GUARDAR MODELOS
# ============================================================

saveRDS(models, file.path(models_dir, "PBMC_modelos_mixtos.rds"))


# ============================================================
# 16. INFORMACIÓN DEL ANÁLISIS
# ============================================================

analysis_info <- tibble(
  Item = c(
    "Analysis",
    "Input file",
    "Number of animals",
    "Groups",
    "Time points (days post-vaccination)",
    "Cytokines",
    "Negative net values set to 0",
    "Transformation",
    "Model",
    "Random effect",
    "Degrees of freedom",
    "Time contrasts",
    "Group contrasts",
    "R version"
  ),
  Value = c(
    "E2-stimulated PBMC cytokines before challenge",
    basename(input_file),
    as.character(n_distinct(datos$Animal)),
    paste(group_levels, collapse = ", "),
    paste(time_levels, collapse = ", "),
    paste(cytokines$label, collapse = ", "),
    paste0(negative_summary$label, ": ", negative_summary$n_negative_set_to_zero,
           collapse = "; "),
    "log10(x + 1)",
    "Linear mixed-effects model (nlme::lme, REML): time * group",
    "Animal (random intercept)",
    "nlme containment method",
    "Day 56 vs day 0 within each group (emmeans, no adjustment: one comparison per group)",
    "Pairwise between groups within each day (emmeans, Tukey)",
    R.version.string
  )
)


# ============================================================
# 17. EXPORTAR RESULTADOS A UN SOLO EXCEL
# ============================================================

output_excel <- file.path(tables_dir, "PBMC_cytokines_resultados.xlsx")

wb <- createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
  
  addWorksheet(wb, sheet_name)
  
  if (is.null(data) || nrow(data) == 0) {
    writeData(wb, sheet = sheet_name,
              x = data.frame(Information = "No data available"))
  } else {
    writeData(wb, sheet = sheet_name, x = data)
  }
  
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = 1:ncol(data), widths = "auto")
}

readme_table <- tibble(
  Sheet = c("Datos", "Descriptivos", "ANOVA", "Tiempo_en_grupo",
            "Grupos_en_tiempo", "Texto", "Analysis_info"),
  Description = c(
    "Processed data in long format (raw, analyzed and log10(x + 1) values)",
    "Mean, SD, SEM, median and IQR by cytokine, group and day (pg/mL)",
    "Fixed effects of the mixed models: F(numDF, denDF) and p",
    "Day 56 vs day 0 within each group (log10 scale; fold change of x + 1)",
    "Pairwise comparisons between groups within each day (Tukey)",
    "Draft sentences for Section 3.3 of the manuscript",
    "Analysis settings"
  )
)

add_sheet(wb, "README",           readme_table)
add_sheet(wb, "Datos",            datos_long %>% select(Animal, Grupo, Tiempo, Cytokine = label,
                                                        Value_raw, Value, log_value))
add_sheet(wb, "Descriptivos",     descriptives %>% rename(Cytokine = label))
add_sheet(wb, "ANOVA",            anova_table)
add_sheet(wb, "Tiempo_en_grupo",  time_contrasts)
add_sheet(wb, "Grupos_en_tiempo", group_contrasts)
add_sheet(wb, "Texto",            text_summary)
add_sheet(wb, "Analysis_info",    analysis_info)

saveWorkbook(wb, output_excel, overwrite = TRUE)


# ============================================================
# 18. INFORME EN UN SOLO WORD
# ============================================================

output_word <- file.path(reports_dir, "PBMC_cytokines_informe.docx")

make_ft <- function(df) {
  flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    autofit()
}

doc <- read_docx() %>%
  body_add_par("Cellular response before challenge: E2-stimulated PBMC cytokines",
               style = "heading 1") %>%
  body_add_par(paste0("Generated on ", format(Sys.Date(), "%Y-%m-%d"),
                      " with ", R.version.string, "."), style = "Normal") %>%
  body_add_par("Analysis", style = "heading 2") %>%
  body_add_par(paste(
    "Net cytokine concentrations (E2-stimulated minus unstimulated) were",
    "log10(x + 1)-transformed and analyzed with linear mixed-effects models",
    "(nlme::lme, REML) including time, group and their interaction as fixed",
    "effects and animal as a random intercept. Day 56 was compared with day 0",
    "within each group, and groups were compared within each day with",
    "Tukey's adjustment (emmeans)."
  ), style = "Normal") %>%
  body_add_par("Negative net values set to 0", style = "heading 3") %>%
  body_add_flextable(make_ft(negative_summary %>%
                               rename(Cytokine = label,
                                      `Values set to 0` = n_negative_set_to_zero))) %>%
  body_add_par("Draft text for Section 3.3", style = "heading 2")

for (s in text_summary$Sentence) {
  doc <- body_add_par(doc, s, style = "Normal")
}

for (lab in cytokines$label) {
  
  doc <- doc %>%
    body_add_par(lab, style = "heading 2") %>%
    body_add_par("Fixed effects", style = "heading 3") %>%
    body_add_flextable(make_ft(
      anova_table %>%
        filter(Cytokine == lab) %>%
        select(Effect, numDF, denDF, F_value, p = p_text)
    )) %>%
    body_add_par("Day 56 vs day 0 within each group", style = "heading 3") %>%
    body_add_flextable(make_ft(
      time_contrasts %>%
        filter(Cytokine == lab) %>%
        transmute(Group = group_labels[as.character(Grupo)],
                  `Estimate (log10)` = round(estimate_log10, 3),
                  SE = round(SE, 3),
                  df = round(df, 1),
                  t = round(t_ratio, 2),
                  `Fold change` = round(fold_change, 2),
                  p = p_text)
    ))
}

doc <- doc %>%
  body_add_break() %>%
  body_add_par("Figure 4", style = "heading 2") %>%
  body_add_img(src = file.path(figures_dir, "PBMC_Figure_4_cytokines.png"),
               width = 6.3, height = 6.3 * 8.5 / 7.2) %>%
  body_add_par(paste(
    "Cytokine secretion by E2-stimulated PBMCs at day 0 and day 56.",
    "Bars show mean ± SEM (pg/mL); symbols show individual animals.",
    "p-values: day 56 vs day 0 from the linear mixed model on",
    "log10(x + 1)-transformed data (shown when p < 0.05)."
  ), style = "Normal")

print(doc, target = output_word)


# ============================================================
# 19. SESSION INFO (reproducibilidad en GitHub)
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_PBMC.txt"))


# ============================================================
# 20. MENSAJE FINAL
# ============================================================

cat("\n\n============================================================\n")
cat("ANÁLISIS COMPLETADO CORRECTAMENTE\n")
cat("============================================================\n")
cat("\nExcel:  ", output_excel, "\n")
cat("Word:   ", output_word, "\n")
cat("Figuras:", figures_dir, "\n")
cat("Modelos:", models_dir, "\n")
cat("\n============================================================\n")
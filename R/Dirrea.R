# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# DIARREA DESPUÉS DEL DESAFÍO
#
# Sección 3.6 (Diarrhea), Figura 8C–D y filas de diarrea de la Tabla 4
#
# Puntaje diario de diarrea (Tabla 2): 0 normal, 1 pastosa/blanda,
# 2 acuosa, 3 hemorrágica; intermedios 0–1 y 1–2 codificados 0.5 y 1.5
# Días post-desafío: 0 a 28
#
# Análisis (Secciones 2.10 y 2.11 del manuscrito):
#   - AUC del puntaje por animal (trapecio, días 0–28)
#   - Kruskal–Wallis + Dunn con ajuste de Benjamini–Hochberg
#   - Sensibilidad: AUC excluyendo los días 0–3
#   - Días-animal con diarrea (puntaje > 0) por grupo
#   - Tabla 4: animales con puntaje > 0 y con puntaje ≥ 2
#
# Salidas (carpeta results/ del proyecto):
#   figures/Diarrhea_Figure_8C_course.(pdf|tiff|png)
#   figures/Diarrhea_Figure_8D_AUC.(pdf|tiff|png)
#   figures/Clinical_Figure_8_complete.(pdf|tiff|png)  (si existe la
#     salida del script de temperatura con los paneles 8A y 8B)
#   tables/Diarrea_resultados.xlsx
#   reports/Diarrea_informe.docx
#   sessionInfo_diarrea.txt
#
# Autor: Santiago Salazar
# ============================================================


# ============================================================
# 1. OPCIONES Y PAQUETES
# ============================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "ragg", "systemfonts", "ggplot2", "dplyr", "tidyr", "tibble",
  "readxl", "openxlsx", "rstatix", "pracma", "patchwork", "cowplot",
  "purrr", "officer", "flextable", "ggpubr"
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
  library(rstatix)
  library(pracma)
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

input_file <- file.path(data_dir, "Diarrea.xlsx")

# Códigos de la planilla para los puntajes intermedios (Sección 2.10):
# "1.2" = entre 1 y 2 -> 1.5 ; "0.1" = entre 0 y 1 -> 0.5
intermediate_codes <- c("1.2" = 1.5, "0.1" = 0.5)

sensitivity_exclude_days <- 0:3   # heces blandas en los días 1–2 (Sección 3.6)

group_levels <- c("Control", "25ug", "50ug", "100ug", "Commercial")

# Paneles 8A y 8B guardados por el script de temperatura (opcional)
temperature_panels_file <- file.path(models_dir, "Temperature_Figure_8AB_plots.rds")


# ============================================================
# 6. LECTURA DE DATOS
# ============================================================

if (!file.exists(input_file)) {
  stop(paste0("No se encontró el archivo:\n", input_file))
}

cat("\nArchivo de entrada:\n", input_file, "\n")

datos_raw <- readxl::read_excel(input_file, sheet = 1)
names(datos_raw) <- trimws(names(datos_raw))

required_columns <- c("DIIO", "Vacuna", "Tiempo", "Score")
missing_columns  <- setdiff(required_columns, names(datos_raw))

if (length(missing_columns) > 0) {
  stop(paste0("Faltan las siguientes columnas:\n",
              paste(missing_columns, collapse = ", ")))
}


# ============================================================
# 7. PREPARACIÓN DE DATOS
# ============================================================

recode_group <- function(x) {
  
  x  <- trimws(as.character(x))
  xl <- tolower(x)
  
  xn <- ifelse(
    grepl("[0-9]+\\s*(ug|µg|mcg)", xl, perl = TRUE),
    sub("^.*?([0-9]+)\\s*(ug|µg|mcg).*$", "\\1", xl, perl = TRUE),
    gsub("[^0-9]", "", xl)
  )
  
  dplyr::case_when(
    toupper(x) == "G1"                         ~ "Control",
    toupper(x) == "G2"                         ~ "25ug",
    toupper(x) == "G3"                         ~ "50ug",
    toupper(x) == "G4"                         ~ "100ug",
    toupper(x) == "G5"                         ~ "Commercial",
    grepl("control|placebo", xl)               ~ "Control",
    grepl("comer|commer|comm|cattle", xl)      ~ "Commercial",
    xn == "25"                                 ~ "25ug",
    xn == "50"                                 ~ "50ug",
    xn == "100"                                ~ "100ug",
    TRUE                                       ~ NA_character_
  )
}

datos <- datos_raw %>%
  mutate(
    Score_raw_text = trimws(gsub(",", ".", as.character(Score))),
    Score_raw      = suppressWarnings(as.numeric(Score_raw_text)),
    recoded        = Score_raw_text %in% names(intermediate_codes),
    Score          = ifelse(recoded, intermediate_codes[Score_raw_text], Score_raw),
    Animal         = factor(as.character(DIIO)),
    Grupo          = factor(recode_group(Vacuna), levels = group_levels),
    Tiempo         = as.numeric(Tiempo)
  ) %>%
  filter(!is.na(Tiempo)) %>%
  arrange(Grupo, Animal, Tiempo)

if (any(is.na(datos$Grupo))) {
  stop("Hay grupos que no se pudieron reconocer: ",
       paste(unique(datos$Vacuna[is.na(datos$Grupo)]), collapse = ", "))
}

invalid_scores <- datos %>% filter(is.na(Score) | !Score %in% c(0, 0.5, 1, 1.5, 2, 2.5, 3))

if (nrow(invalid_scores) > 0) {
  warning("Hay puntajes fuera de la escala 0–3 o vacíos; revisar la hoja 'Puntajes_revisar'.")
}

recode_summary <- datos %>%
  filter(recoded) %>%
  count(Score_raw_text, Score, name = "n_readings")


# ============================================================
# 8. COMPROBACIÓN BÁSICA
# ============================================================

cat("\n================ DATOS =================\n")
cat("\nAnimales por grupo:\n")
print(datos %>% distinct(Animal, Grupo) %>% count(Grupo))
cat("\nDías observados por animal (rango):\n")
print(datos %>% count(Animal) %>% summarise(min = min(n), max = max(n)))
cat("\nPuntajes intermedios recodificados:\n")
print(recode_summary)


# ============================================================
# 9. FUENTE, PALETA, FORMAS Y TEMA (iguales al resto de los scripts)
# ============================================================

base_family <- "Arial"

if (!base_family %in% systemfonts::system_fonts()$family) {
  warning("No se encontró la fuente ", base_family, "; se usará 'sans'.")
  base_family <- "sans"
}

update_geom_defaults("text",  list(family = base_family))
update_geom_defaults("label", list(family = base_family))

palette_nature <- c(
  "Control"    = "#4D4D4D",
  "25ug"       = "#004B87",
  "50ug"       = "#E69F00",
  "100ug"      = "#2ECC71",
  "Commercial" = "#D55E00"
)

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
  
  cowplot::theme_cowplot(font_size = base_size, font_family = base_family) +
    theme(
      plot.title           = element_text(face = "plain", size = base_size, hjust = 0.5),
      axis.title           = element_text(size = base_size),
      axis.text            = element_text(size = base_size, color = "grey15"),
      legend.title         = element_blank(),
      legend.position      = "top",
      legend.justification = "center",
      panel.grid.major.y   = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank()
    )
}


# ============================================================
# 10. FUNCIÓN PARA GUARDAR FIGURAS
# ============================================================

save_nature <- function(plot, filename, width = 7, height = 4) {
  
  # PDF con cairo_pdf y TIFF/PNG con ragg: usan la fuente Arial instalada
  # en el sistema y admiten caracteres Unicode (µ, ±)
  ggsave(file.path(figures_dir, paste0(filename, ".pdf")),
         plot = plot, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf)
  
  ggsave(file.path(figures_dir, paste0(filename, ".tiff")),
         plot = plot, width = width, height = height, units = "in",
         dpi = 600, compression = "lzw", bg = "white",
         device = ragg::agg_tiff)
  
  ggsave(file.path(figures_dir, paste0(filename, ".png")),
         plot = plot, width = width, height = height, units = "in",
         dpi = 600, bg = "white",
         device = ragg::agg_png)
}


# ============================================================
# 11. FORMATO DE VALORES P
# ============================================================

fmt_p_table <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
fmt_p_text  <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))


# ============================================================
# 12. CURSO DIARIO (MEDIA ± SEM)
# ============================================================

daily_summary <- datos %>%
  group_by(Grupo, Tiempo) %>%
  summarise(
    n    = sum(!is.na(Score)),
    mean = mean(Score, na.rm = TRUE),
    sd   = sd(Score, na.rm = TRUE),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

max_after_day3 <- daily_summary %>%
  filter(Tiempo > 3) %>%
  group_by(Grupo) %>%
  slice_max(mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Grupo, day_of_max = Tiempo, max_mean_score = mean)


# ============================================================
# 13. AUC POR ANIMAL (DÍAS 0–28 Y SENSIBILIDAD SIN DÍAS 0–3)
# ============================================================

auc_by_animal <- function(df) {
  df %>%
    filter(!is.na(Score)) %>%
    arrange(Tiempo) %>%
    group_by(Grupo, Animal) %>%
    summarise(
      n_days = n(),
      AUC    = pracma::trapz(Tiempo, Score),
      .groups = "drop"
    )
}

auc_main <- auc_by_animal(datos)
auc_sens <- auc_by_animal(datos %>% filter(!Tiempo %in% sensitivity_exclude_days))

describe_auc <- function(df) {
  df %>%
    group_by(Grupo) %>%
    summarise(
      n      = n(),
      median = median(AUC),
      q1     = quantile(AUC, 0.25),
      q3     = quantile(AUC, 0.75),
      mean   = mean(AUC),
      sd     = sd(AUC),
      .groups = "drop"
    )
}

auc_descriptives      <- describe_auc(auc_main)
auc_descriptives_sens <- describe_auc(auc_sens)


# ============================================================
# 14. KRUSKAL–WALLIS + DUNN (BH)
# ============================================================

kw_main <- rstatix::kruskal_test(auc_main, AUC ~ Grupo)
kw_sens <- rstatix::kruskal_test(auc_sens, AUC ~ Grupo)

dunn_main <- rstatix::dunn_test(auc_main, AUC ~ Grupo, p.adjust.method = "BH")
dunn_sens <- rstatix::dunn_test(auc_sens, AUC ~ Grupo, p.adjust.method = "BH")

kw_table <- bind_rows(
  kw_main %>% mutate(Analysis = "Days 0–28"),
  kw_sens %>% mutate(Analysis = "Excluding days 0–3")
) %>%
  transmute(Analysis, n, H = round(statistic, 2), df, p_value = p, p_text = fmt_p_table(p))

dunn_table <- bind_rows(
  dunn_main %>% mutate(Analysis = "Days 0–28"),
  dunn_sens %>% mutate(Analysis = "Excluding days 0–3")
) %>%
  transmute(Analysis, group1, group2, n1, n2, z = round(statistic, 2),
            p_unadjusted = p, p_BH = p.adj, p_BH_text = fmt_p_table(p.adj))


# ============================================================
# 15. DÍAS-ANIMAL CON DIARREA Y TABLA 4
# ============================================================

animal_days <- datos %>%
  filter(!is.na(Score)) %>%
  group_by(Grupo) %>%
  summarise(
    n_animals       = n_distinct(Animal),
    animal_days     = n(),
    days_with_score = sum(Score > 0),
    pct             = 100 * days_with_score / animal_days,
    .groups = "drop"
  )

table4_diarrhea <- datos %>%
  filter(!is.na(Score)) %>%
  group_by(Grupo, Animal) %>%
  summarise(any_diarrhea = any(Score > 0), any_watery = any(Score >= 2),
            .groups = "drop") %>%
  group_by(Grupo) %>%
  summarise(n = n(),
            `Diarrhea, score > 0`  = sum(any_diarrhea),
            `Diarrhea, score ≥ 2`  = sum(any_watery),
            .groups = "drop")


# ============================================================
# 16. FRASES LISTAS PARA EL MANUSCRITO (Sección 3.6)
# ============================================================

dunn_line <- function(dunn, g1, g2) {
  r <- dunn %>% filter((group1 == g1 & group2 == g2) | (group1 == g2 & group2 == g1))
  r$p.adj
}

desc_line <- auc_descriptives %>%
  mutate(txt = sprintf("%s %.2f (%.2f–%.2f)", group_labels[as.character(Grupo)], median, q1, q3)) %>%
  pull(txt) %>% paste(collapse = "; ")

mean_line <- auc_descriptives %>%
  mutate(txt = sprintf("%s %.2f ± %.2f", group_labels[as.character(Grupo)], mean, sd)) %>%
  pull(txt) %>% paste(collapse = "; ")

other_p <- dunn_main %>%
  filter(!(group1 == "25ug" & group2 %in% c("Control", "Commercial")),
         !(group2 == "25ug" & group1 %in% c("Control", "Commercial")))

days_line <- animal_days %>%
  mutate(txt = sprintf("%s %d of %d (%.0f%%)", group_labels[as.character(Grupo)],
                       days_with_score, animal_days, pct)) %>%
  pull(txt) %>% paste(collapse = "; ")

text_summary <- tibble(
  Item = c("Kruskal–Wallis", "Median (IQR)", "Mean ± SD",
           "25 µg vs commercial", "25 µg vs control (descriptive)",
           "Other comparisons", "Sensitivity (excluding days 0–3)",
           "Animal-days with diarrhea"),
  Sentence = c(
    sprintf("H = %.2f, df = %d, %s.", kw_main$statistic, kw_main$df, fmt_p_text(kw_main$p)),
    paste0(desc_line, "."),
    paste0(mean_line, "."),
    paste0("Dunn-BH ", fmt_p_text(dunn_line(dunn_main, "25ug", "Commercial")), "."),
    paste0("Dunn-BH ", fmt_p_text(dunn_line(dunn_main, "25ug", "Control")), "."),
    paste0("All other adjusted p ≥ ", sprintf("%.3f", min(other_p$p.adj)), "."),
    sprintf("Kruskal–Wallis %s; 25 µg vs commercial, Dunn-BH %s; 25 µg vs control, Dunn-BH %s.",
            fmt_p_text(kw_sens$p),
            fmt_p_text(dunn_line(dunn_sens, "25ug", "Commercial")),
            fmt_p_text(dunn_line(dunn_sens, "25ug", "Control"))),
    paste0(days_line, ".")
  )
)

cat("\n================ RESUMEN PARA EL TEXTO =================\n")
for (i in seq_len(nrow(text_summary))) {
  cat("\n", text_summary$Item[i], ": ", text_summary$Sentence[i], sep = "")
}
cat("\n\nTabla 4 (diarrea):\n")
print(table4_diarrhea)


# ============================================================
# 17. FIGURA 8C: CURSO DIARIO DEL PUNTAJE (MEDIA ± SEM)
# ============================================================

fig8C <- ggplot(daily_summary,
                aes(x = Tiempo, y = mean, color = Grupo, fill = Grupo, group = Grupo)) +
  geom_ribbon(aes(ymin = pmax(0, mean - sem), ymax = mean + sem),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = Grupo), size = 1.6) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_y_continuous(limits = c(0, 3.05), breaks = 0:3,
                     labels = c("0: Normal", "1: Pasty", "2: Watery", "3: Bloody")) +
  scale_x_continuous(breaks = seq(0, 28, 2)) +
  labs(x = "Days after challenge", y = "Diarrhea score (mean ± SEM)") +
  theme_nature(10)

save_nature(fig8C, "Diarrhea_Figure_8C_course", width = 6.5, height = 3.6)


# ============================================================
# 18. FIGURA 8D: AUC (BARRAS SOLO p < 0.05 ENTRE VACUNADOS;
#     EL CONTROL ES DESCRIPTIVO)
# ============================================================

brackets_8D <- dunn_main %>%
  filter(p.adj < 0.05, group1 != "Control", group2 != "Control")

fig8D <- ggplot(auc_main, aes(x = Grupo, y = AUC)) +
  geom_boxplot(aes(color = Grupo, fill = Grupo), width = 0.5,
               outlier.shape = NA, alpha = 0.15, linewidth = 0.4) +
  geom_point(aes(color = Grupo, fill = Grupo, shape = Grupo),
             position = position_jitter(width = 0.1, height = 0, seed = 1),
             size = 1.6, alpha = 0.9) +
  scale_x_discrete(labels = group_labels) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Cumulative diarrhea score (AUC)") +
  guides(color = "none", fill = "none", shape = "none") +
  theme_nature(10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

if (nrow(brackets_8D) > 0) {
  brackets_8D <- brackets_8D %>%
    rstatix::add_y_position(data = auc_main, formula = AUC ~ Grupo, step.increase = 0.1) %>%
    mutate(label = fmt_p_text(p.adj))
  
  fig8D <- fig8D +
    ggpubr::stat_pvalue_manual(brackets_8D, label = "label", tip.length = 0.01,
                               bracket.size = 0.35, size = 2.4)
}

save_nature(fig8D, "Diarrhea_Figure_8D_AUC", width = 2.6, height = 3.8)


# ============================================================
# 19. FIGURA 8 COMPLETA (A–D), SI EXISTEN LOS PANELES DE TEMPERATURA
# ============================================================

if (file.exists(temperature_panels_file)) {
  
  panels_AB <- readRDS(temperature_panels_file)
  
  fig8_complete <- (panels_AB$fig8A + panels_AB$fig8B +
                      plot_layout(widths = c(2.6, 1))) /
    (fig8C + guides(color = "none", fill = "none", shape = "none") + fig8D +
       plot_layout(widths = c(2.6, 1))) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))
  
  save_nature(fig8_complete, "Clinical_Figure_8_complete", width = 9.0, height = 7.4)
  
} else {
  message("No se encontraron los paneles 8A–B; se guardan solo 8C y 8D.\n",
          "Corre primero el script de temperatura para obtener la Figura 8 completa.")
}


# ============================================================
# 20. INFORMACIÓN DEL ANÁLISIS
# ============================================================

analysis_info <- tibble(
  Item = c("Input file", "Number of animals", "Days", "Intermediate scores",
           "AUC", "Tests", "Sensitivity analysis", "R version"),
  Value = c(
    basename(input_file),
    as.character(n_distinct(datos$Animal)),
    paste(range(datos$Tiempo), collapse = "–"),
    paste0(names(intermediate_codes), " -> ", intermediate_codes, collapse = "; "),
    "Trapezoidal AUC of the daily diarrhea score per animal",
    "Kruskal–Wallis and Dunn's test with Benjamini–Hochberg adjustment",
    paste0("AUC excluding days ", paste(range(sensitivity_exclude_days), collapse = "–")),
    R.version.string
  )
)


# ============================================================
# 21. EXPORTAR RESULTADOS A UN SOLO EXCEL
# ============================================================

output_excel <- file.path(tables_dir, "Diarrea_resultados.xlsx")

wb <- createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  if (is.null(data) || nrow(data) == 0) data <- data.frame(Information = "No data available")
  writeData(wb, sheet = sheet_name, x = data)
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(data)), widths = "auto")
}

add_sheet(wb, "Datos",             datos %>% select(Animal, Grupo, Tiempo, Score_raw_text, Score))
add_sheet(wb, "Puntajes_revisar",  invalid_scores %>% select(Animal, Grupo, Tiempo, Score_raw_text))
add_sheet(wb, "Recodificados",     recode_summary)
add_sheet(wb, "Curso_diario",      daily_summary)
add_sheet(wb, "Max_despues_dia3",  max_after_day3)
add_sheet(wb, "AUC_animal",        auc_main)
add_sheet(wb, "AUC_descriptivos",  auc_descriptives)
add_sheet(wb, "AUC_sin_dias0_3",   auc_descriptives_sens)
add_sheet(wb, "Kruskal_Wallis",    kw_table)
add_sheet(wb, "Dunn_BH",           dunn_table)
add_sheet(wb, "Dias_animal",       animal_days)
add_sheet(wb, "Tabla4_diarrea",    table4_diarrhea)
add_sheet(wb, "Texto",             text_summary)
add_sheet(wb, "Analysis_info",     analysis_info)

saveWorkbook(wb, output_excel, overwrite = TRUE)


# ============================================================
# 22. INFORME EN UN SOLO WORD
# ============================================================

output_word <- file.path(reports_dir, "Diarrea_informe.docx")

make_ft <- function(df) {
  flextable(df) %>% theme_booktabs() %>% fontsize(size = 9, part = "all") %>% autofit()
}

doc <- read_docx() %>%
  body_add_par("Diarrhea after challenge", style = "heading 1") %>%
  body_add_par(paste0("Generated on ", format(Sys.Date(), "%Y-%m-%d"),
                      " with ", R.version.string, "."), style = "Normal") %>%
  body_add_par("Draft text for Section 3.6", style = "heading 2")

for (i in seq_len(nrow(text_summary))) {
  doc <- body_add_par(doc, paste0(text_summary$Item[i], ": ", text_summary$Sentence[i]),
                      style = "Normal")
}

doc <- doc %>%
  body_add_par("AUC by group", style = "heading 2") %>%
  body_add_flextable(make_ft(auc_descriptives %>%
                               transmute(Group = group_labels[as.character(Grupo)], n,
                                         `Median (IQR)` = sprintf("%.2f (%.2f–%.2f)", median, q1, q3),
                                         `Mean ± SD` = sprintf("%.2f ± %.2f", mean, sd)))) %>%
  body_add_par("Kruskal–Wallis", style = "heading 3") %>%
  body_add_flextable(make_ft(kw_table %>% select(Analysis, n, H, df, p = p_text))) %>%
  body_add_par("Dunn's test (BH)", style = "heading 3") %>%
  body_add_flextable(make_ft(dunn_table %>%
                               mutate(group1 = group_labels[group1], group2 = group_labels[group2]) %>%
                               select(Analysis, group1, group2, z, p = p_BH_text))) %>%
  body_add_par("Table 4: diarrhea rows", style = "heading 2") %>%
  body_add_flextable(make_ft(table4_diarrhea %>%
                               mutate(Grupo = group_labels[as.character(Grupo)]) %>%
                               rename(Group = Grupo))) %>%
  body_add_break() %>%
  body_add_par("Figure 8C–D", style = "heading 2") %>%
  body_add_img(file.path(figures_dir, "Diarrhea_Figure_8C_course.png"),
               width = 6.3, height = 6.3 * 3.6 / 6.5) %>%
  body_add_img(file.path(figures_dir, "Diarrhea_Figure_8D_AUC.png"),
               width = 2.4, height = 2.4 * 3.8 / 2.6)

if (file.exists(file.path(figures_dir, "Clinical_Figure_8_complete.png"))) {
  doc <- doc %>%
    body_add_par("Figure 8 (complete)", style = "heading 2") %>%
    body_add_img(file.path(figures_dir, "Clinical_Figure_8_complete.png"),
                 width = 6.5, height = 6.5 * 7.4 / 9.0)
}

print(doc, target = output_word)


# ============================================================
# 23. SESSION INFO Y MENSAJE FINAL
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_diarrea.txt"))

cat("\n\n============================================================\n")
cat("ANÁLISIS COMPLETADO CORRECTAMENTE\n")
cat("============================================================\n")
cat("\nExcel:  ", output_excel, "\n")
cat("Word:   ", output_word, "\n")
cat("Figuras:", figures_dir, "\n")
cat("\n============================================================\n")
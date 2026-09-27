# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# CITOQUINAS SÉRICAS DESPUÉS DEL DESAFÍO
#
# Sección 3.5, Figura 6 y Figura 7 del manuscrito
#
# Citoquinas: IFN-γ, IL-17, TNF-α, IL-6 e IL-4
# Días post-desafío: 0, 3, 6, 14 y 21
#
# Análisis (Secciones 2.9 y 2.11 del manuscrito):
#   - Transformación log10(x + 1)
#   - Modelo lineal mixto (lmerTest::lmer): Grupo * Tiempo,
#     intercepto aleatorio por animal; F con gl de Satterthwaite
#   - Grupos dentro de cada día (emmeans, ajuste de Tukey)
#   - Cada día vs día 0 dentro de cada grupo (emmeans, Dunnett; exploratorio)
#   - AUC (pg/mL × días) y pico por animal: Kruskal–Wallis + Dunn (BH)
#   - log2 fold change individual vs día 0 (Figura 7, descriptivo)
#
# Salidas (carpeta results/ del proyecto):
#   figures/Serum_Figure_6_cytokine_kinetics.pdf / .tiff / .png
#   figures/Serum_Figure_7_cytokine_log2FC_heatmap.pdf / .tiff / .png
#   figures/Serum_Figure_S_cytokine_AUC.pdf / .tiff / .png
#   tables/Serum_cytokines_resultados.xlsx
#   reports/Serum_cytokines_informe.docx
#   models/Serum_cytokines_modelos_mixtos.rds
#   sessionInfo_serum_cytokines.txt
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
  "lme4",
  "lmerTest",
  "emmeans",
  "rstatix",
  "pracma",
  "patchwork",
  "cowplot",
  "purrr",
  "scales",
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
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(rstatix)
  library(pracma)
  library(patchwork)
  library(cowplot)
  library(purrr)
  library(scales)
  library(officer)
  library(flextable)
})

# Grados de libertad de Satterthwaite también en emmeans (coherente con Métodos)
emm_options(lmer.df = "satterthwaite")


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

input_file  <- file.path(data_dir, "Citoquinas_desafio.xlsx")
input_sheet <- "Desafio"

# Columna del Excel y etiqueta (orden de los paneles de la Figura 6)
cytokines <- tibble::tribble(
  ~var,          ~label,
  "IFN_pg_ml",   "IFN-γ",
  "IL17_pg_ml",  "IL-17",
  "TNF_pg_mL",   "TNF-α",
  "IL6_pg_mL",   "IL-6",
  "IL4_pg_mL",   "IL-4"
)

# Orden de filas en el heatmap (de arriba hacia abajo, como en la figura actual)
heatmap_row_order <- c("IL-4", "IL-6", "TNF-α", "IL-17", "IFN-γ")

# Constante que se suma antes de calcular el fold change individual
# (evita dividir por 0 cuando la concentración basal es 0).
# Con basales cercanos a 0, el log2FC puede ser muy grande: interpretar con cautela.
fc_pseudocount <- 0.01

alpha_level  <- 0.05
group_levels <- c("Control", "25ug", "50ug", "100ug", "Commercial")
time_levels  <- c(0, 3, 6, 14, 21)


# ============================================================
# 6. LECTURA DE DATOS
# ============================================================

if (!file.exists(input_file)) {
  stop(paste0("No se encontró el archivo:\n", input_file))
}

cat("\nArchivo de entrada:\n", input_file, "\n")

datos_raw <- readxl::read_excel(input_file, sheet = input_sheet)

required_columns <- c("Animal", "Grupo", "Tiempo", cytokines$var)
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
  
  # Dosis: número que precede a "ug"/"µg" (evita tomar el "2" de "rE2 25ug")
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
    Animal  = factor(Animal),
    Grupo   = factor(recode_group(Grupo), levels = group_levels),
    Tiempo  = as.numeric(Tiempo),
    TiempoF = factor(Tiempo, levels = time_levels)
  )

if (any(is.na(datos$Grupo))) {
  stop("Hay grupos que no se pudieron reconocer: ",
       paste(unique(datos_raw$Grupo[is.na(datos$Grupo)]), collapse = ", "))
}

if (any(is.na(datos$TiempoF))) {
  stop("Hay días distintos de 0, 3, 6, 14 y 21: ",
       paste(unique(datos$Tiempo[is.na(datos$TiempoF)]), collapse = ", "))
}

datos_long <- datos %>%
  select(Animal, Grupo, Tiempo, TiempoF, all_of(cytokines$var)) %>%
  pivot_longer(cols = all_of(cytokines$var), names_to = "var", values_to = "Value") %>%
  mutate(Value = as.numeric(Value)) %>%
  left_join(cytokines, by = "var") %>%
  mutate(
    Cytokine  = factor(label, levels = cytokines$label),
    log_value = log10(Value + 1)
  ) %>%
  select(-label) %>%
  arrange(Cytokine, Grupo, Animal, Tiempo)

if (any(datos_long$Value < 0, na.rm = TRUE)) {
  stop("Hay concentraciones negativas en el archivo; revisar los datos.")
}


# ============================================================
# 8. COMPROBACIÓN BÁSICA
# ============================================================

cat("\n================ DATOS =================\n")

cat("\nAnimales por grupo:\n")
print(datos %>% distinct(Animal, Grupo) %>% count(Grupo))

cat("\nObservaciones por grupo y día:\n")
print(datos %>% count(Grupo, Tiempo) %>%
        pivot_wider(names_from = Tiempo, values_from = n))

cat("\nValores faltantes por citoquina:\n")
print(datos_long %>% group_by(Cytokine) %>% summarise(missing = sum(is.na(Value))))


# ============================================================
# 9. PALETA, FORMAS Y TEMA (iguales al resto de los scripts)
# ============================================================

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
  
  cowplot::theme_cowplot(font_size = base_size, font_family = "sans") +
    theme(
      plot.title          = element_text(face = "plain", size = base_size, hjust = 0.5),
      axis.title          = element_text(size = base_size),
      axis.text           = element_text(size = base_size, color = "grey15"),
      strip.background    = element_rect(fill = "grey95", color = NA),
      strip.text          = element_text(face = "bold", size = base_size),
      legend.title        = element_blank(),
      legend.position     = "top",
      legend.justification = "center",
      panel.grid.major.y  = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank()
    )
}


# ============================================================
# 10. FUNCIÓN PARA GUARDAR FIGURAS
# ============================================================

save_nature <- function(plot, filename, width = 7, height = 4) {
  
  # cairo_pdf: incrusta las fuentes y admite caracteres Unicode (γ, α, µ, ×)
  ggsave(file.path(figures_dir, paste0(filename, ".pdf")),
         plot = plot, width = width, height = height, units = "in",
         device = grDevices::cairo_pdf)
  
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

fmt_p_table <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
fmt_p_text  <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))


# ============================================================
# 12. ESTADÍSTICA DESCRIPTIVA (pg/mL)
# ============================================================

descriptives <- datos_long %>%
  group_by(Cytokine, Grupo, Tiempo) %>%
  summarise(
    n         = sum(!is.na(Value)),
    mean      = mean(Value, na.rm = TRUE),
    sd        = sd(Value, na.rm = TRUE),
    sem       = sd / sqrt(n),
    median    = median(Value, na.rm = TRUE),
    q1        = quantile(Value, 0.25, na.rm = TRUE),
    q3        = quantile(Value, 0.75, na.rm = TRUE),
    geom_mean = 10^mean(log_value, na.rm = TRUE) - 1,
    .groups   = "drop"
  )


# ============================================================
# 13. MODELOS LINEALES MIXTOS POR CITOQUINA
# ============================================================

cytokine_data <- split(datos_long, datos_long$Cytokine)[cytokines$label]

models <- purrr::map(cytokine_data, function(df) {
  lmerTest::lmer(log_value ~ Grupo * TiempoF + (1 | Animal), data = df, REML = TRUE)
})


# ------------------------------------------------------------
# 13.1 Efectos fijos: F(NumDF, DenDF) y p (tipo III, Satterthwaite)
# ------------------------------------------------------------

anova_table <- purrr::imap_dfr(models, function(m, lab) {
  
  a <- as.data.frame(anova(m, type = 3, ddf = "Satterthwaite"))
  
  tibble(
    Cytokine = lab,
    Effect   = rownames(a),
    NumDF    = a$NumDF,
    DenDF    = round(a$DenDF, 1),
    F_value  = round(a$`F value`, 2),
    p_value  = a$`Pr(>F)`,
    p_text   = fmt_p_table(a$`Pr(>F)`)
  ) %>%
    mutate(Effect = recode(Effect,
                           "Grupo"         = "Group",
                           "TiempoF"       = "Time",
                           "Grupo:TiempoF" = "Group × time"))
})


# ------------------------------------------------------------
# 13.2 Grupos dentro de cada día (Tukey)
# ------------------------------------------------------------

group_contrasts <- purrr::imap_dfr(models, function(m, lab) {
  
  as.data.frame(pairs(emmeans(m, ~ Grupo | TiempoF), adjust = "tukey")) %>%
    as_tibble() %>%
    transmute(
      Cytokine       = lab,
      Day            = as.numeric(as.character(TiempoF)),
      Contrast       = contrast,
      estimate_log10 = estimate,
      SE, df,
      t_ratio        = t.ratio,
      p_value        = p.value,
      p_text         = fmt_p_table(p.value),
      significant    = p.value < alpha_level
    )
})


# ------------------------------------------------------------
# 13.3 Cada día vs día 0 dentro de cada grupo (Dunnett; exploratorio)
# ------------------------------------------------------------

time_contrasts <- purrr::imap_dfr(models, function(m, lab) {
  
  as.data.frame(contrast(emmeans(m, ~ TiempoF | Grupo),
                         method = "trt.vs.ctrl", ref = 1)) %>%
    as_tibble() %>%
    transmute(
      Cytokine       = lab,
      Grupo,
      Contrast       = contrast,
      estimate_log10 = estimate,
      SE, df,
      t_ratio        = t.ratio,
      p_value        = p.value,
      p_text         = fmt_p_table(p.value)
    )
})


# ============================================================
# 14. MÉTRICAS POR ANIMAL: AUC Y PICO
# ============================================================

animal_metrics <- datos_long %>%
  arrange(Tiempo) %>%
  group_by(Cytokine, Grupo, Animal) %>%
  summarise(
    n_points     = sum(!is.na(Value)),
    baseline     = Value[Tiempo == 0][1],
    peak         = max(Value, na.rm = TRUE),
    time_to_peak = Tiempo[which.max(Value)][1],
    auc          = if (n_points >= 2) {
      pracma::trapz(Tiempo[!is.na(Value)], Value[!is.na(Value)])
    } else {
      NA_real_
    },
    .groups = "drop"
  )

metric_summary <- animal_metrics %>%
  group_by(Cytokine, Grupo) %>%
  summarise(
    n         = n(),
    auc_mean  = mean(auc, na.rm = TRUE),
    auc_sd    = sd(auc, na.rm = TRUE),
    auc_median = median(auc, na.rm = TRUE),
    peak_mean = mean(peak, na.rm = TRUE),
    peak_sd   = sd(peak, na.rm = TRUE),
    peak_median = median(peak, na.rm = TRUE),
    .groups   = "drop"
  )

# Kruskal–Wallis y Dunn (BH), como indica la Sección 2.11 para AUC y pico de TNF-α
metric_long <- animal_metrics %>%
  select(Cytokine, Grupo, Animal, auc, peak) %>%
  pivot_longer(c(auc, peak), names_to = "Metric", values_to = "Value") %>%
  filter(!is.na(Value))

kw_table <- metric_long %>%
  group_by(Cytokine, Metric) %>%
  rstatix::kruskal_test(Value ~ Grupo) %>%
  ungroup() %>%
  transmute(Cytokine, Metric, n, H = round(statistic, 2), df,
            p_value = p, p_text = fmt_p_table(p))

dunn_table <- metric_long %>%
  group_by(Cytokine, Metric) %>%
  rstatix::dunn_test(Value ~ Grupo, p.adjust.method = "BH") %>%
  ungroup() %>%
  transmute(Cytokine, Metric, group1, group2, n1, n2,
            z = round(statistic, 2), p_unadjusted = p, p_BH = p.adj,
            p_BH_text = fmt_p_table(p.adj))


# ============================================================
# 15. LOG2 FOLD CHANGE INDIVIDUAL VS DÍA 0 (FIGURA 7)
# ============================================================

log2fc_individual <- datos_long %>%
  filter(Tiempo > 0) %>%
  left_join(
    datos_long %>%
      filter(Tiempo == 0) %>%
      select(Animal, Cytokine, baseline = Value),
    by = c("Animal", "Cytokine")
  ) %>%
  mutate(log2_fc = log2((Value + fc_pseudocount) / (baseline + fc_pseudocount)))

# Media de los log2FC individuales = log2 de la media geométrica del fold change
heatmap_data <- log2fc_individual %>%
  group_by(Cytokine, Grupo, Tiempo) %>%
  summarise(
    n       = sum(!is.na(log2_fc)),
    log2_fc = mean(log2_fc, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 16. FRASES LISTAS PARA EL MANUSCRITO (Sección 3.5)
# ============================================================

fmt_F <- function(row) {
  paste0("F", row$NumDF, ",", round(row$DenDF), " = ",
         sprintf("%.2f", row$F_value), ", ", fmt_p_text(row$p_value))
}

text_summary <- purrr::map_dfr(cytokines$label, function(lab) {
  
  a <- anova_table %>% filter(Cytokine == lab)
  
  tibble(
    Cytokine = lab,
    Sentence = paste0(
      lab, ": time, ",          fmt_F(a %>% filter(Effect == "Time")),
      "; group, ",              fmt_F(a %>% filter(Effect == "Group")),
      "; group × time, ",       fmt_F(a %>% filter(Effect == "Group × time")), "."
    )
  )
})

# Comparaciones significativas entre grupos, por día
significant_group_text <- group_contrasts %>%
  filter(significant) %>%
  mutate(Sentence = paste0(Cytokine, ", day ", Day, ": ", Contrast, ", ",
                           fmt_p_text(p_value), ".")) %>%
  select(Cytokine, Sentence)

cat("\n================ RESUMEN PARA EL TEXTO =================\n")
for (s in text_summary$Sentence) cat("\n", s)
cat("\n\nComparaciones significativas entre grupos (Tukey):\n")
for (s in significant_group_text$Sentence) cat("\n", s)
cat("\n")


# ============================================================
# 17. FIGURA 6: CINÉTICA DE CITOQUINAS SÉRICAS
#     Puntos: animales; líneas: media; barras: ± SEM
# ============================================================

facet_labels <- setNames(paste0(cytokines$label, " (pg/mL)"), cytokines$label)

fig6_kinetics <- ggplot() +
  
  geom_point(
    data = datos_long,
    aes(x = Tiempo, y = Value, color = Grupo, fill = Grupo, shape = Grupo),
    position = position_jitter(width = 0.35, height = 0, seed = 1),
    size = 1.0, alpha = 0.35
  ) +
  
  geom_line(
    data = descriptives,
    aes(x = Tiempo, y = mean, color = Grupo, group = Grupo),
    linewidth = 0.7
  ) +
  
  geom_errorbar(
    data = descriptives,
    aes(x = Tiempo, ymin = pmax(0, mean - sem), ymax = mean + sem, color = Grupo),
    width = 0.5, linewidth = 0.35
  ) +
  
  facet_wrap(~ Cytokine, scales = "free_y", ncol = 3,
             labeller = as_labeller(facet_labels)) +
  
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_x_continuous(breaks = time_levels) +
  
  labs(x = "Days after challenge", y = "Concentration (pg/mL)") +
  
  theme_nature(10) +
  theme(axis.text.x = element_text(size = 8))

save_nature(fig6_kinetics, "Serum_Figure_6_cytokine_kinetics", width = 9.0, height = 6.5)


# ============================================================
# 18. FIGURA 7: HEATMAP DE LOG2 FOLD CHANGE (EN INGLÉS)
# ============================================================

fig7_heatmap <- heatmap_data %>%
  mutate(Cytokine = factor(Cytokine, levels = rev(heatmap_row_order))) %>%
  ggplot(aes(x = factor(Tiempo), y = Cytokine, fill = log2_fc)) +
  
  geom_tile(color = "white", linewidth = 0.8) +
  
  geom_text(
    aes(label = sub("^-0\\.00$", "0.00", sprintf("%.2f", log2_fc)),
        color = ifelse(abs(log2_fc) > 3, "white", "black")),
    size = 3.0, fontface = "bold"
  ) +
  scale_color_identity() +
  
  facet_wrap(~ Grupo, nrow = 1, labeller = as_labeller(group_labels)) +
  
  scale_fill_gradient2(
    low = "#002147", mid = "#FFFFFF", high = "#680018",
    midpoint = 0, limits = c(-5, 5), oob = scales::squish,
    name = "log2 fold change\n(vs. day 0)"
  ) +
  
  labs(x = "Days after challenge", y = NULL) +
  
  theme_cowplot(font_size = 10) +
  theme(
    axis.text        = element_text(color = "black"),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(0.3, "lines"),
    legend.title     = element_text(face = "bold", size = 9)
  )

save_nature(fig7_heatmap, "Serum_Figure_7_cytokine_log2FC_heatmap", width = 10.5, height = 4.5)


# ============================================================
# 19. FIGURA SUPLEMENTARIA: AUC POR GRUPO
# ============================================================

figS_auc <- ggplot(animal_metrics, aes(x = Grupo, y = auc)) +
  geom_boxplot(aes(fill = Grupo), width = 0.45, outlier.shape = NA,
               alpha = 0.2, linewidth = 0.4) +
  geom_point(aes(color = Grupo, fill = Grupo, shape = Grupo),
             position = position_jitter(width = 0.12, height = 0, seed = 1),
             size = 1.5, alpha = 0.85) +
  facet_wrap(~ Cytokine, scales = "free_y", ncol = 3,
             labeller = as_labeller(setNames(cytokines$label, cytokines$label))) +
  scale_x_discrete(labels = group_labels) +
  scale_color_manual(values = palette_nature, labels = group_labels) +
  scale_fill_manual(values = palette_nature, labels = group_labels) +
  scale_shape_manual(values = shape_groups, labels = group_labels) +
  labs(x = NULL, y = "AUC (pg/mL × days)") +
  theme_nature(10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

save_nature(figS_auc, "Serum_Figure_S_cytokine_AUC", width = 8.5, height = 5.5)


# ============================================================
# 20. GUARDAR MODELOS
# ============================================================

saveRDS(models, file.path(models_dir, "Serum_cytokines_modelos_mixtos.rds"))


# ============================================================
# 21. INFORMACIÓN DEL ANÁLISIS
# ============================================================

analysis_info <- tibble(
  Item = c(
    "Analysis", "Input file", "Number of animals", "Groups",
    "Days after challenge", "Cytokines", "Transformation", "Model",
    "Random effect", "Degrees of freedom", "Group contrasts",
    "Time contrasts", "AUC and peak", "log2 fold change", "R version"
  ),
  Value = c(
    "Serum cytokines after challenge",
    paste0(basename(input_file), " (sheet: ", input_sheet, ")"),
    as.character(n_distinct(datos$Animal)),
    paste(group_levels, collapse = ", "),
    paste(time_levels, collapse = ", "),
    paste(cytokines$label, collapse = ", "),
    "log10(x + 1)",
    "Linear mixed-effects model (lmerTest::lmer, REML): group * time",
    "Animal (random intercept)",
    "Satterthwaite (type III F tests and emmeans)",
    "Pairwise between groups within each day (emmeans, Tukey)",
    "Each day vs day 0 within each group (emmeans, Dunnett; exploratory)",
    "Trapezoidal AUC (pg/mL × days) and peak per animal; Kruskal–Wallis and Dunn's test (BH)",
    paste0("Mean of individual log2[(x + ", fc_pseudocount, ") / (day 0 + ",
           fc_pseudocount, ")]"),
    R.version.string
  )
)


# ============================================================
# 22. EXPORTAR RESULTADOS A UN SOLO EXCEL
# ============================================================

output_excel <- file.path(tables_dir, "Serum_cytokines_resultados.xlsx")

wb <- createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
  
  addWorksheet(wb, sheet_name)
  
  if (is.null(data) || nrow(data) == 0) {
    data <- data.frame(Information = "No data available")
  }
  
  writeData(wb, sheet = sheet_name, x = data)
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(data)), widths = "auto")
}

readme_table <- tibble(
  Sheet = c("Datos", "Descriptivos", "ANOVA", "Grupos_en_dia", "Dias_vs_dia0",
            "Metricas_animal", "Metricas_resumen", "Kruskal_Wallis", "Dunn_BH",
            "Log2FC_heatmap", "Log2FC_individual", "Texto", "Analysis_info"),
  Description = c(
    "Processed data in long format",
    "Mean, SD, SEM, median, IQR and geometric mean by cytokine, group and day (pg/mL)",
    "Fixed effects of the mixed models: F(NumDF, DenDF) and p",
    "Pairwise comparisons between groups within each day (Tukey)",
    "Each day vs day 0 within each group (Dunnett; exploratory)",
    "AUC (pg/mL × days), peak and time to peak per animal",
    "AUC and peak summary by group",
    "Kruskal–Wallis tests for AUC and peak",
    "Dunn's pairwise tests with BH adjustment for AUC and peak",
    "Mean log2 fold change vs day 0 (values shown in Figure 7)",
    "Individual log2 fold change vs day 0",
    "Draft sentences for Section 3.5 of the manuscript",
    "Analysis settings"
  )
)

add_sheet(wb, "README",            readme_table)
add_sheet(wb, "Datos",             datos_long %>% select(Animal, Grupo, Tiempo, Cytokine, Value, log_value))
add_sheet(wb, "Descriptivos",      descriptives)
add_sheet(wb, "ANOVA",             anova_table)
add_sheet(wb, "Grupos_en_dia",     group_contrasts)
add_sheet(wb, "Dias_vs_dia0",      time_contrasts)
add_sheet(wb, "Metricas_animal",   animal_metrics)
add_sheet(wb, "Metricas_resumen",  metric_summary)
add_sheet(wb, "Kruskal_Wallis",    kw_table)
add_sheet(wb, "Dunn_BH",           dunn_table)
add_sheet(wb, "Log2FC_heatmap",    heatmap_data)
add_sheet(wb, "Log2FC_individual", log2fc_individual %>%
            select(Animal, Grupo, Tiempo, Cytokine, Value, baseline, log2_fc))
add_sheet(wb, "Texto",             bind_rows(text_summary, significant_group_text))
add_sheet(wb, "Analysis_info",     analysis_info)

saveWorkbook(wb, output_excel, overwrite = TRUE)


# ============================================================
# 23. INFORME EN UN SOLO WORD
# ============================================================

output_word <- file.path(reports_dir, "Serum_cytokines_informe.docx")

make_ft <- function(df) {
  flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    autofit()
}

doc <- read_docx() %>%
  body_add_par("Serum cytokines after challenge", style = "heading 1") %>%
  body_add_par(paste0("Generated on ", format(Sys.Date(), "%Y-%m-%d"),
                      " with ", R.version.string, "."), style = "Normal") %>%
  body_add_par("Analysis", style = "heading 2") %>%
  body_add_par(paste(
    "Serum cytokine concentrations were log10(x + 1)-transformed and analyzed",
    "with linear mixed-effects models (lmerTest::lmer, REML) including group,",
    "time and their interaction as fixed effects and animal as a random",
    "intercept. Fixed effects were tested with type III F tests and",
    "Satterthwaite's degrees of freedom; groups were compared within each day",
    "with Tukey's adjustment (emmeans). Cumulative exposure (trapezoidal AUC)",
    "and peak concentration per animal were compared with the Kruskal–Wallis",
    "test followed by Dunn's test with Benjamini–Hochberg adjustment."
  ), style = "Normal") %>%
  body_add_par("Draft text for Section 3.5", style = "heading 2")

for (s in c(text_summary$Sentence, significant_group_text$Sentence)) {
  doc <- body_add_par(doc, s, style = "Normal")
}

doc <- doc %>%
  body_add_par("Fixed effects of the mixed models", style = "heading 2") %>%
  body_add_flextable(make_ft(
    anova_table %>% select(Cytokine, Effect, NumDF, DenDF, F_value, p = p_text)
  )) %>%
  body_add_par("Significant between-group comparisons (Tukey)", style = "heading 2") %>%
  body_add_flextable(make_ft(
    group_contrasts %>%
      filter(significant) %>%
      transmute(Cytokine, Day, Contrast,
                `Estimate (log10)` = round(estimate_log10, 3),
                t = round(t_ratio, 2), p = p_text)
  )) %>%
  body_add_par("AUC and peak: Kruskal–Wallis", style = "heading 2") %>%
  body_add_flextable(make_ft(
    kw_table %>% select(Cytokine, Metric, n, H, df, p = p_text)
  )) %>%
  body_add_par("AUC and peak by group (mean ± SD)", style = "heading 2") %>%
  body_add_flextable(make_ft(
    metric_summary %>%
      transmute(Cytokine, Group = group_labels[as.character(Grupo)], n,
                AUC = sprintf("%.2f ± %.2f", auc_mean, auc_sd),
                Peak = sprintf("%.2f ± %.2f", peak_mean, peak_sd))
  )) %>%
  body_add_break() %>%
  body_add_par("Figure 6", style = "heading 2") %>%
  body_add_img(src = file.path(figures_dir, "Serum_Figure_6_cytokine_kinetics.png"),
               width = 6.5, height = 6.5 * 6.5 / 9.0) %>%
  body_add_par(paste(
    "Kinetics of serum cytokines after challenge. Symbols show individual",
    "animals, lines show group means and error bars show SEM. Each panel has",
    "an independent y-axis."
  ), style = "Normal") %>%
  body_add_par("Figure 7", style = "heading 2") %>%
  body_add_img(src = file.path(figures_dir, "Serum_Figure_7_cytokine_log2FC_heatmap.png"),
               width = 6.5, height = 6.5 * 4.5 / 10.5) %>%
  body_add_par(paste(
    "Relative change in serum cytokine concentration at 3, 6, 14 and 21 days",
    "after challenge. Values are the mean of individual log2 fold changes",
    "relative to day 0, equivalent to the log2 of the geometric mean fold change."
  ), style = "Normal")

print(doc, target = output_word)


# ============================================================
# 24. SESSION INFO (reproducibilidad en GitHub)
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_serum_cytokines.txt"))


# ============================================================
# 25. MENSAJE FINAL
# ============================================================

cat("\n\n============================================================\n")
cat("ANÁLISIS COMPLETADO CORRECTAMENTE\n")
cat("============================================================\n")
cat("\nExcel:  ", output_excel, "\n")
cat("Word:   ", output_word, "\n")
cat("Figuras:", figures_dir, "\n")
cat("Modelos:", models_dir, "\n")
cat("\n============================================================\n")
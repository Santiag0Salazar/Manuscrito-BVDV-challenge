# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# SIGNOS CLÍNICOS Y PUNTAJE CLÍNICO COMPUESTO DESPUÉS DEL DESAFÍO
#
# Sección 3.6 ("Other clinical signs and composite score"), Tabla 4,
# AUC del puntaje de temperatura y figura suplementaria del puntaje
# compuesto
#
# Entrada: data/processed/Signos_clinicos.xlsx, generado por
# 00_Limpieza_datos_clinicos.R a partir de la planilla de registro diario
# (una hoja por día, "13-mayo" ... "10-junio"). Puntajes de la Tabla 2 ya
# codificados: intermedios 0-1 = 0.5 y 1-2 = 1.5; anotaciones solo
# descriptivas = 0 (listadas en la hoja "Anotaciones_texto").
#
# Análisis (Secciones 2.10 y 2.11):
#   - AUC (trapecio, días 0–28) de cada signo, del puntaje compuesto
#     (suma diaria de todos los signos excepto temperatura) y del
#     puntaje de temperatura
#   - Kruskal–Wallis + Dunn con ajuste de Benjamini–Hochberg
#   - Contribución de cada signo al AUC del puntaje compuesto
#   - Tabla 4: número de animales con cada signo (al menos una vez)
#   - Fisher exacto: recombinantes (agrupados) vs comercial
#
# Salidas (carpeta results/ del proyecto):
#   figures/Composite_score_Figure_S2.(pdf|tiff|png)
#   tables/Signos_clinicos_resultados.xlsx
#   tables/Clinical_scores_AUC_per_animal.csv  (para las correlaciones, 3.7)
#   reports/Signos_clinicos_informe.docx
#   sessionInfo_signos_clinicos.txt
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

# Tabla limpia generada por 00_Limpieza_datos_clinicos.R (puntajes ya
# codificados: intermedios 0.5/1.5 y anotaciones de texto = 0)
input_file <- file.path(data_dir, "Signos_clinicos.xlsx")

challenge_date <- as.Date("2026-05-13")   # día 0 (desafío)
study_days     <- 0:28
pm_days        <- 0:2                     # lectura vespertina de temperatura (Sección 2.10)
peak_day       <- 8                       # día de temperatura máxima (Tabla 4)

# Número de la figura suplementaria del puntaje compuesto (Figura S2)
composite_fig_id   <- "S2"
composite_fig_file <- paste0("Composite_score_Figure_", composite_fig_id)

# Signos de la Tabla 2 (columnas de Signos_clinicos.xlsx)
sign_columns <- tribble(
  ~pattern,                 ~Sign,
  "^desc\\.?\\s*nasal",     "Nasal discharge",
  "^desc\\.?\\s*ocular",    "Ocular discharge",
  "^tos",                   "Cough",
  "^disnea",                "Dyspnea",
  "^comport",               "Behavior",
  "^diarrea",               "Diarrhea",
  "^apetito",               "Appetite",
  "^dh$",                   "Dehydration",
  "^t.*rectal",             "Temperature"
)

composite_signs <- setdiff(sign_columns$Sign, "Temperature")

group_levels <- c("Control", "25ug", "50ug", "100ug", "Commercial")


# ============================================================
# 6. LECTURA DE DATOS (TABLA LIMPIA)
# ============================================================

if (!file.exists(input_file)) {
  stop(paste0("No se encontró el archivo:\n", input_file,
              "\nCorre primero 00_Limpieza_datos_clinicos.R"))
}

cat("\nArchivo de entrada:\n", input_file, "\n")

daily_raw        <- readxl::read_excel(input_file, sheet = "Puntajes")
text_log         <- readxl::read_excel(input_file, sheet = "Anotaciones_texto")
intermediate_log <- readxl::read_excel(input_file, sheet = "Intermedios")

required_columns <- c("Animal", "Grupo", "Day", sign_columns$Sign, "Temperature_pm")
missing_columns  <- setdiff(required_columns, names(daily_raw))
if (length(missing_columns) > 0) {
  stop("Faltan columnas en ", basename(input_file), ": ", paste(missing_columns, collapse = ", "))
}

text_log <- text_log %>%
  mutate(Animal = factor(as.character(Animal)),
         Grupo  = factor(Grupo, levels = group_levels))


# ============================================================
# 7. COMPROBACIÓN DE PUNTAJES
# ============================================================

allowed_scores <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)

scores_check <- daily_raw %>%
  select(Animal, Grupo, Day, all_of(sign_columns$Sign)) %>%
  pivot_longer(all_of(sign_columns$Sign), names_to = "Sign", values_to = "score") %>%
  filter(is.na(score) | !score %in% allowed_scores |
           (Sign == "Temperature" & !score %in% c(0, 2, 3)))

if (nrow(scores_check) > 0) {
  warning("Hay puntajes faltantes o fuera de escala; revisar la hoja 'Puntajes_revisar'.")
}


# ============================================================
# 8. TABLA DIARIA POR ANIMAL, PUNTAJE COMPUESTO Y TEMPERATURA
# ============================================================

daily <- daily_raw %>%
  mutate(
    Animal = factor(as.character(Animal)),
    Grupo  = factor(Grupo, levels = group_levels),
    Day    = as.numeric(Day),
    across(all_of(c(sign_columns$Sign, "Temperature_pm")), as.numeric),
    # Puntaje diario de temperatura: máximo de la mañana y la tarde en
    # los días con lectura vespertina (días 0–2)
    Temperature_pm    = ifelse(Day %in% pm_days, Temperature_pm, NA_real_),
    Temperature_daily = pmax(Temperature, Temperature_pm, na.rm = TRUE),
    Composite         = rowSums(across(all_of(composite_signs)), na.rm = TRUE)
  ) %>%
  arrange(Grupo, Animal, Day)

if (any(is.na(daily$Grupo))) stop("Grupos no reconocidos en ", basename(input_file))

cat("\n================ DATOS =================\n")
cat("\nAnimales por grupo:\n")
print(daily %>% distinct(Animal, Grupo) %>% count(Grupo))
cat("\nDías observados por animal (rango):\n")
print(daily %>% count(Animal) %>% summarise(min = min(n), max = max(n)))
cat("\nAnotaciones de texto contadas como ausentes:",
    sum(text_log$type == "descriptive (absent)"), "\n")


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

palette_signs <- c(
  "Diarrhea"         = "#5E3C99",
  "Cough"            = "#B2ABD2",
  "Nasal discharge"  = "#FDB863",
  "Ocular discharge" = "#E66101",
  "Dyspnea"          = "#80CDC1",
  "Behavior"         = "#018571",
  "Appetite"         = "#A6611A",
  "Dehydration"      = "#BABABA"
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
# 12. AUC POR ANIMAL (DÍAS 0–28)
# ============================================================

auc_vars <- c(composite_signs, "Composite", "Temperature_daily", "Temperature")

auc_animal <- daily %>%
  arrange(Animal, Day) %>%
  group_by(Grupo, Animal) %>%
  summarise(
    n_days = n(),
    across(all_of(auc_vars), ~ pracma::trapz(Day, .x)),
    .groups = "drop"
  ) %>%
  rename(Temperature_AUC_am_pm = Temperature_daily,
         Temperature_AUC_am    = Temperature)

auc_long <- auc_animal %>%
  pivot_longer(-c(Grupo, Animal, n_days), names_to = "Variable", values_to = "AUC")

describe_auc <- function(df) {
  df %>%
    group_by(Variable, Grupo) %>%
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

auc_descriptives <- describe_auc(auc_long)


# ============================================================
# 13. KRUSKAL–WALLIS + DUNN (BH)
#     Puntaje compuesto (primario) y puntaje de temperatura
# ============================================================

tested_vars <- c(
  "Composite"             = "Composite score (excluding temperature)",
  "Diarrhea"              = "Diarrhea score (cross-check with Diarrea_post_desafio.R)",
  "Temperature_AUC_am_pm" = "Temperature score (daily maximum of morning and evening readings)",
  "Temperature_AUC_am"    = "Temperature score (morning readings only)"
)

kw_table <- purrr::map_dfr(names(tested_vars), function(v) {
  auc_long %>%
    filter(Variable == v) %>%
    rstatix::kruskal_test(AUC ~ Grupo) %>%
    mutate(Variable = v, Description = tested_vars[[v]])
}) %>%
  transmute(Variable, Description, n, H = round(statistic, 2), df,
            p_value = p, p_text = fmt_p_table(p))

dunn_all <- purrr::map_dfr(names(tested_vars), function(v) {
  auc_long %>%
    filter(Variable == v) %>%
    rstatix::dunn_test(AUC ~ Grupo, p.adjust.method = "BH") %>%
    mutate(Variable = v)
})

dunn_table <- dunn_all %>%
  transmute(Variable, group1, group2, n1, n2, z = round(statistic, 2),
            p_unadjusted = p, p_BH = p.adj, p_BH_text = fmt_p_table(p.adj))

kw_row   <- function(v) kw_table %>% filter(Variable == v)
dunn_p   <- function(v, g1, g2) {
  dunn_all %>%
    filter(Variable == v,
           (group1 == g1 & group2 == g2) | (group1 == g2 & group2 == g1)) %>%
    pull(p.adj)
}


# ============================================================
# 14. CONTRIBUCIÓN DE CADA SIGNO AL PUNTAJE COMPUESTO
# ============================================================

sign_contribution <- auc_long %>%
  filter(Variable %in% composite_signs) %>%
  group_by(Grupo, Sign = Variable) %>%
  summarise(total_AUC = sum(AUC), mean_AUC = mean(AUC), .groups = "drop") %>%
  group_by(Grupo) %>%
  mutate(pct_of_composite = 100 * total_AUC / sum(total_AUC)) %>%
  ungroup() %>%
  mutate(Sign = factor(Sign, levels = names(palette_signs)))

diarrhea_share <- sign_contribution %>%
  filter(Sign == "Diarrhea") %>%
  select(Grupo, pct_of_composite)


# ============================================================
# 15. TABLA 4: ANIMALES CON CADA SIGNO (AL MENOS UNA VEZ)
#     Temperatura: puntaje de la mañana (nota de la Tabla 4)
# ============================================================

per_animal <- daily %>%
  group_by(Grupo, Animal) %>%
  summarise(
    `Diarrhea, score > 0`          = any(Diarrhea > 0),
    `Diarrhea, score ≥ 2`          = any(Diarrhea >= 2),
    `Temperature ≥ 39.2 °C`        = any(Temperature >= 2),
    `Temperature > 40 °C`          = any(Temperature >= 3),
    `Temperature > 40 °C on day 8` = any(Temperature[Day == peak_day] >= 3),
    `Nasal discharge`              = any(`Nasal discharge` > 0),
    `Ocular discharge`             = any(`Ocular discharge` > 0),
    Cough                          = any(Cough > 0),
    Dyspnea                        = any(Dyspnea > 0),
    `Altered behavior`        = any(Behavior > 0),
    `Reduced appetite`             = any(Appetite > 0),
    Dehydration                    = any(Dehydration > 0),
    .groups = "drop"
  )

names(per_animal) <- sub("on day 8", paste("on day", peak_day), names(per_animal))

table4 <- per_animal %>%
  group_by(Grupo) %>%
  summarise(across(-Animal, ~ sum(.x)), n = n(), .groups = "drop") %>%
  pivot_longer(-c(Grupo, n), names_to = "Sign", values_to = "Animals") %>%
  mutate(Group = paste0(group_labels[as.character(Grupo)], " (n = ", n, ")")) %>%
  select(Sign, Group, Animals) %>%
  pivot_wider(names_from = Group, values_from = Animals) %>%
  mutate(Sign = factor(Sign, levels = setdiff(names(per_animal), c("Grupo", "Animal")))) %>%
  arrange(Sign) %>%
  mutate(Sign = as.character(Sign))


# ============================================================
# 16. FISHER EXACTO: RECOMBINANTES (AGRUPADOS) VS COMERCIAL
# ============================================================

fisher_rec_vs_comm <- function(variable) {

  d <- per_animal %>%
    filter(Grupo != "Control") %>%
    mutate(Set = ifelse(Grupo == "Commercial", "Commercial", "Recombinant"),
           y   = .data[[variable]])

  tab <- table(factor(d$Set, levels = c("Recombinant", "Commercial")),
               factor(d$y, levels = c(TRUE, FALSE)))

  tibble(
    Outcome     = variable,
    Recombinant = paste0(tab["Recombinant", "TRUE"], "/", sum(tab["Recombinant", ])),
    Commercial  = paste0(tab["Commercial", "TRUE"], "/", sum(tab["Commercial", ])),
    p_value     = fisher.test(tab)$p.value,
    p_text      = fmt_p_table(fisher.test(tab)$p.value)
  )
}

day_peak_var <- paste("Temperature > 40 °C on day", peak_day)

per_animal <- per_animal %>%
  left_join(
    daily %>% filter(Day == peak_day) %>%
      transmute(Animal, fever_peak_day = Temperature >= 2),
    by = "Animal"
  )

fisher_table <- bind_rows(
  fisher_rec_vs_comm("fever_peak_day") %>%
    mutate(Outcome = paste("Temperature ≥ 39.2 °C on day", peak_day)),
  fisher_rec_vs_comm(day_peak_var),
  fisher_rec_vs_comm("Temperature > 40 °C"),
  fisher_rec_vs_comm("Temperature ≥ 39.2 °C")
)


# ============================================================
# 17. DETALLE DE LOS OTROS SIGNOS (DÍAS CON SIGNO POR ANIMAL)
# ============================================================

other_signs <- c("Nasal discharge", "Ocular discharge", "Cough", "Dyspnea",
                 "Behavior", "Appetite", "Dehydration")

sign_days_animal <- daily %>%
  select(Grupo, Animal, Day, all_of(other_signs)) %>%
  pivot_longer(all_of(other_signs), names_to = "Sign", values_to = "score") %>%
  group_by(Sign, Grupo, Animal) %>%
  summarise(
    days_observed = n(),
    days_with_sign = sum(score > 0, na.rm = TRUE),
    max_score      = max(score, na.rm = TRUE),
    days           = paste(Day[score > 0 & !is.na(score)], collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(days_with_sign > 0) %>%
  arrange(Sign, Grupo, desc(days_with_sign))

sign_days_group <- daily %>%
  select(Grupo, Animal, Day, all_of(other_signs)) %>%
  pivot_longer(all_of(other_signs), names_to = "Sign", values_to = "score") %>%
  group_by(Sign, Grupo) %>%
  summarise(
    animals_with_sign = n_distinct(Animal[score > 0 & !is.na(score)]),
    n_animals         = n_distinct(Animal),
    animal_days       = sum(score > 0, na.rm = TRUE),
    total_animal_days = n(),
    max_score         = max(score, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 18. FRASES LISTAS PARA EL MANUSCRITO (Sección 3.6)
# ============================================================

comp_desc <- auc_descriptives %>% filter(Variable == "Composite")

desc_line <- comp_desc %>%
  mutate(txt = sprintf("%s %.2f (%.2f–%.2f)", group_labels[as.character(Grupo)], median, q1, q3)) %>%
  pull(txt) %>% paste(collapse = "; ")

mean_line <- comp_desc %>%
  mutate(txt = sprintf("%s %.2f ± %.2f", group_labels[as.character(Grupo)], mean, sd)) %>%
  pull(txt) %>% paste(collapse = "; ")

share_line <- diarrhea_share %>%
  mutate(txt = sprintf("%s %.0f%%", group_labels[as.character(Grupo)], pct_of_composite)) %>%
  pull(txt) %>% paste(collapse = "; ")

other_p_comp <- dunn_all %>%
  filter(Variable == "Composite",
         !(group1 == "25ug" & group2 %in% c("Control", "Commercial")),
         !(group2 == "25ug" & group1 %in% c("Control", "Commercial")))

groups_with <- function(sign_label) {
  r <- table4 %>% filter(Sign == sign_label) %>% select(-Sign)
  paste(paste0(names(r), ": ", unlist(r[1, ])), collapse = "; ")
}

top_animal <- function(sign_label) {
  r <- sign_days_animal %>% filter(Sign == sign_label) %>% slice_max(days_with_sign, n = 1, with_ties = FALSE)
  if (nrow(r) == 0) return("none")
  sprintf("ID%s (%s) on %d of %d days (days %s)", r$Animal, group_labels[as.character(r$Grupo)],
          r$days_with_sign, r$days_observed, r$days)
}

behavior_detail <- sign_days_animal %>%
  filter(Sign %in% c("Behavior", "Appetite")) %>%
  mutate(txt = sprintf("%s: ID%s (%s), days %s, max score %.1f", Sign, Animal,
                       group_labels[as.character(Grupo)], days, max_score)) %>%
  pull(txt) %>% paste(collapse = "; ")

kw_c  <- kw_row("Composite")
kw_t  <- kw_row("Temperature_AUC_am_pm")
kw_ta <- kw_row("Temperature_AUC_am")

text_summary <- tibble(
  Item = c(
    "Nasal discharge (animals)", "Ocular discharge (animals)", "Cough (animals)",
    "Animal with most days of cough", "Behavior and appetite",
    "Diarrhea share of composite AUC",
    "Composite score, Kruskal–Wallis", "Composite score, median (IQR)",
    "Composite score, mean ± SD",
    "Composite: 25 µg vs control (descriptive)", "Composite: 25 µg vs commercial",
    "Composite: other comparisons",
    "Temperature score AUC (am/pm maximum)", "Temperature score AUC (morning only)",
    "Descriptive entries"
  ),
  Sentence = c(
    groups_with("Nasal discharge"),
    groups_with("Ocular discharge"),
    groups_with("Cough"),
    top_animal("Cough"),
    ifelse(behavior_detail == "", "None recorded", behavior_detail),
    paste0(share_line, "."),
    sprintf("H = %.2f, df = %d, %s.", kw_c$H, kw_c$df, fmt_p_text(kw_c$p_value)),
    paste0(desc_line, "."),
    paste0(mean_line, "."),
    paste0("Dunn-BH ", fmt_p_text(dunn_p("Composite", "25ug", "Control")), "."),
    paste0("Dunn-BH ", fmt_p_text(dunn_p("Composite", "25ug", "Commercial")), "."),
    paste0("All other adjusted p ≥ ", sprintf("%.3f", min(other_p_comp$p.adj)), "."),
    sprintf("Kruskal–Wallis H = %.2f, df = %d, %s.", kw_t$H, kw_t$df, fmt_p_text(kw_t$p_value)),
    sprintf("Kruskal–Wallis H = %.2f, df = %d, %s.", kw_ta$H, kw_ta$df, fmt_p_text(kw_ta$p_value)),
    sprintf("%d text entries counted as absent (see sheet 'Anotaciones_texto').",
            sum(text_log$type == "descriptive (absent)"))
  )
)

cat("\n================ RESUMEN PARA EL TEXTO =================\n")
for (i in seq_len(nrow(text_summary))) {
  cat("\n", text_summary$Item[i], ": ", text_summary$Sentence[i], sep = "")
}
cat("\n\nTabla 4:\n")
print(as.data.frame(table4))
cat("\nFisher (recombinantes vs comercial):\n")
print(as.data.frame(fisher_table))


# ============================================================
# 19. FIGURA SUPLEMENTARIA DEL PUNTAJE COMPUESTO
#     (A) curso diario, (B) AUC, (C) contribución de cada signo
# ============================================================

composite_daily <- daily %>%
  group_by(Grupo, Day) %>%
  summarise(
    n    = sum(!is.na(Composite)),
    mean = mean(Composite, na.rm = TRUE),
    sd   = sd(Composite, na.rm = TRUE),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

figS_A <- ggplot(composite_daily,
                 aes(x = Day, y = mean, color = Grupo, fill = Grupo, group = Grupo)) +
  geom_ribbon(aes(ymin = pmax(0, mean - sem), ymax = mean + sem),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = Grupo), size = 1.6) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 28, 2)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Days after challenge", y = "Composite clinical score\n(mean ± SEM)") +
  theme_nature(10)

auc_comp <- auc_long %>% filter(Variable == "Composite")

brackets_S <- dunn_all %>%
  filter(Variable == "Composite", p.adj < 0.05, group1 != "Control", group2 != "Control")

figS_B <- ggplot(auc_comp, aes(x = Grupo, y = AUC)) +
  geom_boxplot(aes(color = Grupo, fill = Grupo), width = 0.5,
               outlier.shape = NA, alpha = 0.15, linewidth = 0.4) +
  geom_point(aes(color = Grupo, fill = Grupo, shape = Grupo),
             position = position_jitter(width = 0.1, height = 0, seed = 1),
             size = 1.6, alpha = 0.9) +
  scale_x_discrete(labels = group_labels) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Cumulative composite score (AUC)") +
  guides(color = "none", fill = "none", shape = "none") +
  theme_nature(10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

if (nrow(brackets_S) > 0) {
  brackets_S <- brackets_S %>%
    rstatix::add_y_position(data = auc_comp, formula = AUC ~ Grupo, step.increase = 0.1) %>%
    mutate(label = fmt_p_text(p.adj))

  figS_B <- figS_B +
    ggpubr::stat_pvalue_manual(brackets_S, label = "label", tip.length = 0.01,
                               bracket.size = 0.35, size = 2.4)
}

contribution_plot_data <- sign_contribution %>%
  filter(total_AUC > 0) %>%
  droplevels()

figS_C <- ggplot(contribution_plot_data,
                 aes(x = Grupo, y = pct_of_composite, fill = Sign)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.2) +
  scale_x_discrete(labels = group_labels) +
  scale_fill_manual(values = palette_signs) +
  scale_y_continuous(limits = c(0, 100.01), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Share of composite AUC (%)") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_nature(10) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"),
        axis.text.x = element_text(angle = 40, hjust = 1))

figS_composite <- figS_A /
  (figS_B + figS_C + plot_layout(widths = c(1, 1.25))) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))

save_nature(figS_composite, composite_fig_file, width = 7.0, height = 6.6)


# ============================================================
# 20. INFORMACIÓN DEL ANÁLISIS
# ============================================================

analysis_info <- tibble(
  Item = c("Input file", "Number of animals", "Days",
           "Coding", "Composite score",
           "Temperature score", "AUC", "Tests", "Table 4", "R version"),
  Value = c(
    basename(input_file),
    as.character(n_distinct(daily$Animal)),
    paste(range(daily$Day), collapse = "–"),
    "Scores coded in 00_Limpieza_datos_clinicos.R: intermediate 0-1 = 0.5 and 1-2 = 1.5; text-only entries = 0 (see Anotaciones_texto)",
    paste("Daily sum of:", paste(composite_signs, collapse = ", ")),
    paste0("Daily maximum of morning and evening scores on days ",
           paste(range(pm_days), collapse = "–"), "; morning score on the other days"),
    "Trapezoidal AUC of the daily score per animal (days 0–28)",
    "Kruskal–Wallis and Dunn's test with Benjamini–Hochberg adjustment; Fisher's exact test",
    "Animals with a score > 0 at least once; temperature rows from the morning score",
    R.version.string
  )
)


# ============================================================
# 21. EXPORTAR RESULTADOS (UN SOLO EXCEL + CSV PARA CORRELACIONES)
# ============================================================

output_excel <- file.path(tables_dir, "Signos_clinicos_resultados.xlsx")

wb <- createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  if (is.null(data) || nrow(data) == 0) data <- data.frame(Information = "No data available")
  writeData(wb, sheet = sheet_name, x = data)
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(data)), widths = "auto")
}

add_sheet(wb, "Datos_diarios",       daily)
add_sheet(wb, "Anotaciones_texto",   text_log)
add_sheet(wb, "Intermedios",         intermediate_log)
add_sheet(wb, "Puntajes_revisar",    scores_check)
add_sheet(wb, "AUC_animal",          auc_animal)
add_sheet(wb, "AUC_descriptivos",    auc_descriptives)
add_sheet(wb, "Kruskal_Wallis",      kw_table)
add_sheet(wb, "Dunn_BH",             dunn_table)
add_sheet(wb, "Contribucion_signos", sign_contribution)
add_sheet(wb, "Tabla4",              table4)
add_sheet(wb, "Tabla4_por_animal",   per_animal)
add_sheet(wb, "Fisher",              fisher_table)
add_sheet(wb, "Signos_por_animal",   sign_days_animal)
add_sheet(wb, "Signos_por_grupo",    sign_days_group)
add_sheet(wb, "Curso_compuesto",     composite_daily)
add_sheet(wb, "Texto",               text_summary)
add_sheet(wb, "Analysis_info",       analysis_info)

saveWorkbook(wb, output_excel, overwrite = TRUE)

# AUC por animal para el script de correlaciones (Sección 3.7)
write.csv(auc_animal, file.path(tables_dir, "Clinical_scores_AUC_per_animal.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")


# ============================================================
# 22. INFORME EN UN SOLO WORD
# ============================================================

output_word <- file.path(reports_dir, "Signos_clinicos_informe.docx")

make_ft <- function(df) {
  flextable(df) %>% theme_booktabs() %>% fontsize(size = 9, part = "all") %>% autofit()
}

doc <- read_docx() %>%
  body_add_par("Clinical signs and composite clinical score after challenge", style = "heading 1") %>%
  body_add_par(paste0("Generated on ", format(Sys.Date(), "%Y-%m-%d"),
                      " with ", R.version.string, "."), style = "Normal") %>%
  body_add_par("Draft text for Section 3.6", style = "heading 2")

for (i in seq_len(nrow(text_summary))) {
  doc <- body_add_par(doc, paste0(text_summary$Item[i], ": ", text_summary$Sentence[i]),
                      style = "Normal")
}

doc <- doc %>%
  body_add_par("Table 4", style = "heading 2") %>%
  body_add_flextable(make_ft(table4)) %>%
  body_add_par("Temperature rows use the morning score (Table 2). Text-only entries were counted as absent.",
               style = "Normal") %>%
  body_add_par("Fisher's exact tests (recombinant pooled vs commercial)", style = "heading 3") %>%
  body_add_flextable(make_ft(fisher_table %>% select(Outcome, Recombinant, Commercial, p = p_text))) %>%
  body_add_par("Composite score AUC by group", style = "heading 2") %>%
  body_add_flextable(make_ft(comp_desc %>%
                               transmute(Group = group_labels[as.character(Grupo)], n,
                                         `Median (IQR)` = sprintf("%.2f (%.2f–%.2f)", median, q1, q3),
                                         `Mean ± SD` = sprintf("%.2f ± %.2f", mean, sd)))) %>%
  body_add_par("Kruskal–Wallis", style = "heading 3") %>%
  body_add_flextable(make_ft(kw_table %>% select(Description, n, H, df, p = p_text))) %>%
  body_add_par("Dunn's test (BH)", style = "heading 3") %>%
  body_add_flextable(make_ft(dunn_table %>%
                               filter(Variable %in% c("Composite", "Temperature_AUC_am_pm")) %>%
                               mutate(group1 = group_labels[group1], group2 = group_labels[group2]) %>%
                               select(Variable, group1, group2, z, p = p_BH_text))) %>%
  body_add_par("Share of the composite AUC by sign (%)", style = "heading 3") %>%
  body_add_flextable(make_ft(sign_contribution %>%
                               filter(total_AUC > 0) %>%
                               mutate(Group = group_labels[as.character(Grupo)],
                                      pct = sprintf("%.1f", pct_of_composite)) %>%
                               select(Group, Sign, pct) %>%
                               pivot_wider(names_from = Group, values_from = pct, values_fill = "0.0"))) %>%
  body_add_par("Other clinical signs by animal", style = "heading 3") %>%
  body_add_flextable(make_ft(sign_days_animal %>%
                               mutate(Group = group_labels[as.character(Grupo)]) %>%
                               select(Sign, Group, Animal, days_with_sign, days_observed, max_score, days))) %>%
  body_add_par("Text-only entries (counted as absent unless coded)", style = "heading 3") %>%
  body_add_flextable(make_ft(text_log %>%
                               mutate(Group = group_labels[as.character(Grupo)]) %>%
                               select(Day, Animal, Group, Sign, raw_text, score))) %>%
  body_add_break() %>%
  body_add_par(paste0("Figure ", composite_fig_id, ": composite clinical score"), style = "heading 2") %>%
  body_add_img(file.path(figures_dir, paste0(composite_fig_file, ".png")),
               width = 6.3, height = 6.3 * 6.6 / 7.0)

print(doc, target = output_word)


# ============================================================
# 23. SESSION INFO Y MENSAJE FINAL
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_signos_clinicos.txt"))

cat("\n\n============================================================\n")
cat("ANÁLISIS COMPLETADO CORRECTAMENTE\n")
cat("============================================================\n")
cat("\nExcel:  ", output_excel, "\n")
cat("CSV:    ", file.path(tables_dir, "Clinical_scores_AUC_per_animal.csv"), "\n")
cat("Word:   ", output_word, "\n")
cat("Figuras:", figures_dir, "\n")
cat("\n============================================================\n")

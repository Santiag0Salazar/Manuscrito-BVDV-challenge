# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# VIREMIA DESPUÉS DEL DESAFÍO
#
# Sección 3.4, Tabla 3 y Figura 5 del manuscrito
#
# Días post-desafío: 0, 3, 6, 14 y 21 (qPCR en suero)
#
# Análisis (Secciones 2.8 y 2.11 del manuscrito):
#   - Carga viral: log10(copias/mL + 1); muestras < LOQ = 0 copias/mL
#   - AUC (días 0–21, trapecio), pico y número de días positivos por animal
#   - AUC: ANOVA de una vía + Tukey (primario); Kruskal–Wallis + Dunn (BH)
#     como análisis de sensibilidad
#   - Pico y días positivos: Kruskal–Wallis + Dunn (BH)
#   - Prevalencia: GLMM binomial penalizado (blme::bglmer), grupo, día
#     (categórico) e interacción, animal aleatorio; prior normal (SD 2.5)
#     en efectos fijos y gamma en la covarianza; días 3, 6 y 14.
#     OR por pares promediados sobre los días, IC 95 % y p con ajuste Tukey
#   - Prevalencia por día: Fisher exacto por pares con ajuste BH
#
# Salidas (carpeta results/ del proyecto):
#   figures/Viremia_Figure_5.pdf / .tiff / .png
#   tables/Viremia_resultados.xlsx
#   reports/Viremia_informe.docx
#   models/Viremia_modelos.rds
#   sessionInfo_viremia.txt
#
# Autor: Santiago Salazar
# ============================================================


# ============================================================
# 1. OPCIONES Y PAQUETES
# ============================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "tibble", "readxl", "openxlsx",
  "lme4", "blme", "emmeans", "rstatix", "ggpubr", "pracma",
  "patchwork", "cowplot", "purrr", "scales", "officer", "flextable"
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
  library(blme)
  library(emmeans)
  library(rstatix)
  library(ggpubr)
  library(pracma)
  library(patchwork)
  library(cowplot)
  library(purrr)
  library(scales)
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

input_file <- file.path(data_dir, "Datos_viremia.xlsx")

# Definición de muestra positiva (Sección 2.8):
# "positive when amplification was detected, irrespective of whether the
#  viral load could be quantified". Las muestras < LOQ tienen 0 copias/mL,
# así que 'Copias_mL > 0' NO las cuenta como positivas.
# Si el Excel tiene una columna 0/1 con la amplificación detectada, indicar
# su nombre aquí; si es NULL, se usa Copias_mL > 0 y se avisa en la consola.
positivity_column <- NULL   # por ejemplo: "Amplificacion" o "Positivo"

glmm_days  <- c(3, 6, 14)   # días 0 y 21: todos negativos (Sección 2.11)
fisher_days <- c(3, 6, 14)

prior_sd     <- 2.5
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

datos_raw <- readxl::read_excel(input_file, sheet = 1)

required_columns <- c("Animal", "Grupo", "Tiempo", "Copias_mL", positivity_column)
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
    Animal    = factor(Animal),
    Grupo     = factor(recode_group(Grupo), levels = group_levels),
    Tiempo    = as.numeric(Tiempo),
    TiempoF   = factor(Tiempo, levels = time_levels),
    Copias_mL = as.numeric(Copias_mL),
    log_load  = log10(Copias_mL + 1)
  )

if (is.null(positivity_column)) {
  message("\nAVISO: no se indicó columna de amplificación; se usa Copias_mL > 0.\n",
          "Las muestras positivas bajo el LOQ (0 copias/mL) quedan como negativas.\n")
  datos$Viremic <- as.integer(datos$Copias_mL > 0)
} else {
  datos$Viremic <- as.integer(as.numeric(datos[[positivity_column]]) > 0)
}

if (any(is.na(datos$Grupo))) {
  stop("Hay grupos que no se pudieron reconocer: ",
       paste(unique(datos_raw$Grupo[is.na(datos$Grupo)]), collapse = ", "))
}

if (any(is.na(datos$TiempoF))) {
  stop("Hay días distintos de 0, 3, 6, 14 y 21: ",
       paste(unique(datos$Tiempo[is.na(datos$TiempoF)]), collapse = ", "))
}

datos <- datos %>% arrange(Grupo, Animal, Tiempo)


# ============================================================
# 8. COMPROBACIÓN BÁSICA
# ============================================================

cat("\n================ DATOS =================\n")

cat("\nAnimales por grupo:\n")
print(datos %>% distinct(Animal, Grupo) %>% count(Grupo))

cat("\nAnimales virémicos por grupo y día:\n")
print(datos %>%
        group_by(Grupo, Tiempo) %>%
        summarise(pos = paste0(sum(Viremic), "/", n()), .groups = "drop") %>%
        pivot_wider(names_from = Tiempo, values_from = pos))


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
  
  cowplot::theme_cowplot(font_size = base_size, font_family = "arial") +
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

fmt_p_table <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
fmt_p_text  <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))


# ============================================================
# 12. RESUMEN POR GRUPO Y DÍA
# ============================================================

viremia_summary <- datos %>%
  group_by(Grupo, Tiempo) %>%
  summarise(
    n           = n(),
    n_viremic   = sum(Viremic),
    pct_viremic = 100 * mean(Viremic),
    mean_log    = mean(log_load),
    sd_log      = sd(log_load),
    sem_log     = sd_log / sqrt(n),
    mean_copies = mean(Copias_mL),
    sd_copies   = sd(Copias_mL),
    .groups     = "drop"
  )


# ============================================================
# 13. MÉTRICAS POR ANIMAL
# ============================================================

animal_metrics <- datos %>%
  arrange(Tiempo) %>%
  group_by(Grupo, Animal) %>%
  summarise(
    auc_log10      = pracma::trapz(Tiempo, log_load),
    peak_log10     = max(log_load),
    peak_copies    = max(Copias_mL),
    viremic_days   = sum(Viremic[Tiempo > 0]),
    time_to_peak   = if (max(Copias_mL) > 0) Tiempo[which.max(Copias_mL)][1] else NA_real_,
    .groups        = "drop"
  )

metric_summary <- animal_metrics %>%
  group_by(Grupo) %>%
  summarise(
    n                  = n(),
    auc_mean           = mean(auc_log10),
    auc_sd             = sd(auc_log10),
    auc_median         = median(auc_log10),
    peak_log10_mean    = mean(peak_log10),
    peak_log10_sd      = sd(peak_log10),
    peak_copies_mean   = mean(peak_copies),
    peak_copies_sd     = sd(peak_copies),
    viremic_days_mean  = mean(viremic_days),
    viremic_days_sd    = sd(viremic_days),
    .groups            = "drop"
  )


# ============================================================
# 14. AUC: ANOVA + TUKEY (PRIMARIO) Y KRUSKAL–WALLIS (SENSIBILIDAD)
# ============================================================

fit_auc <- lm(auc_log10 ~ Grupo, data = animal_metrics)

auc_anova <- as.data.frame(anova(fit_auc)) %>%
  rownames_to_column("Effect") %>%
  as_tibble()

auc_assumptions <- tibble(
  Test    = c("Shapiro–Wilk (residuals)", "Levene"),
  p_value = c(shapiro.test(residuals(fit_auc))$p.value,
              rstatix::levene_test(animal_metrics, auc_log10 ~ Grupo)$p)
)

auc_tukey <- as.data.frame(pairs(emmeans(fit_auc, ~ Grupo), adjust = "tukey")) %>%
  as_tibble() %>%
  mutate(p_text = fmt_p_table(p.value))

# Tabla para las barras de la figura (mismo Tukey, formato rstatix)
auc_tukey_plot <- animal_metrics %>%
  rstatix::tukey_hsd(auc_log10 ~ Grupo)


# ============================================================
# 15. KRUSKAL–WALLIS + DUNN (BH): AUC, PICO Y DÍAS POSITIVOS
# ============================================================

metric_long <- animal_metrics %>%
  select(Grupo, Animal, auc_log10, peak_log10, viremic_days) %>%
  pivot_longer(c(auc_log10, peak_log10, viremic_days),
               names_to = "Metric", values_to = "Value")

kw_table <- metric_long %>%
  group_by(Metric) %>%
  rstatix::kruskal_test(Value ~ Grupo) %>%
  ungroup() %>%
  transmute(Metric, n, H = round(statistic, 2), df, p_value = p, p_text = fmt_p_table(p))

dunn_table <- metric_long %>%
  group_by(Metric) %>%
  rstatix::dunn_test(Value ~ Grupo, p.adjust.method = "BH") %>%
  ungroup()

dunn_for_plot <- function(metric) {
  animal_metrics %>%
    rstatix::dunn_test(as.formula(paste(metric, "~ Grupo")), p.adjust.method = "BH")
}


# ============================================================
# 16. PREVALENCIA: GLMM BINOMIAL PENALIZADO (TABLA 3)
# ============================================================

datos_glmm <- datos %>%
  filter(Tiempo %in% glmm_days) %>%
  mutate(TiempoF = droplevels(TiempoF))

n_coefs <- ncol(model.matrix(~ Grupo * TiempoF, data = datos_glmm))

fit_glmm <- blme::bglmer(
  Viremic ~ Grupo * TiempoF + (1 | Animal),
  data        = datos_glmm,
  family      = binomial,
  cov.prior   = gamma,
  fixef.prior = normal(cov = diag(prior_sd^2, n_coefs))
)

# OR por pares, promediados sobre los días 3, 6 y 14
emm_glmm   <- emmeans(fit_glmm, ~ Grupo, type = "response")
or_pairs   <- pairs(emm_glmm, adjust = "tukey")

or_table <- as.data.frame(summary(or_pairs)) %>%
  as_tibble() %>%
  left_join(
    as.data.frame(confint(or_pairs, adjust = "tukey")) %>%
      as_tibble() %>%
      select(contrast, asymp.LCL, asymp.UCL),
    by = "contrast"
  ) %>%
  transmute(
    Comparison = gsub("ug", " µg", gsub(" / ", " vs ", contrast)),
    OR         = round(odds.ratio, 3),
    CI_95      = sprintf("%.3f–%.3f", asymp.LCL, asymp.UCL),
    p_adjusted = round(p.value, 3),
    p_text     = fmt_p_table(p.value)
  )


# ============================================================
# 17. PREVALENCIA POR DÍA: FISHER POR PARES (BH DENTRO DE CADA DÍA)
# ============================================================

fisher_by_day <- purrr::map_dfr(fisher_days, function(t) {
  
  df_t  <- datos %>% filter(Tiempo == t)
  combs <- combn(group_levels, 2, simplify = FALSE)
  
  purrr::map_dfr(combs, function(pair) {
    
    sub_df <- df_t %>% filter(Grupo %in% pair) %>% droplevels()
    tab    <- table(factor(sub_df$Grupo, levels = pair),
                    factor(sub_df$Viremic, levels = c(0, 1)))
    
    tibble(
      Day    = t,
      group1 = pair[1],
      group2 = pair[2],
      pos1   = paste0(tab[1, "1"], "/", sum(tab[1, ])),
      pos2   = paste0(tab[2, "1"], "/", sum(tab[2, ])),
      p      = fisher.test(tab)$p.value
    )
  })
}) %>%
  group_by(Day) %>%
  mutate(p_BH = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(p_BH_text = fmt_p_table(p_BH))


# ============================================================
# 18. FRASES LISTAS PARA EL MANUSCRITO (Sección 3.4)
# ============================================================

auc_F <- auc_anova %>% filter(Effect == "Grupo")

vs_comm <- function(tab, col_p) {
  tab %>%
    filter(grepl("Commercial", contrast)) %>%
    mutate(txt = paste0(contrast, ", ", fmt_p_text(.data[[col_p]]))) %>%
    pull(txt) %>%
    paste(collapse = "; ")
}

kw_txt <- function(metric) {
  k <- kw_table %>% filter(Metric == metric)
  paste0("Kruskal–Wallis, H = ", k$H, ", ", fmt_p_text(k$p_value))
}

dunn_txt <- function(metric) {
  dunn_table %>%
    filter(Metric == metric, group1 == "Commercial" | group2 == "Commercial") %>%
    mutate(txt = paste0(group1, " vs ", group2, ", Dunn-BH ", fmt_p_text(p.adj))) %>%
    pull(txt) %>%
    paste(collapse = "; ")
}

min_fisher <- fisher_by_day %>% slice_min(p_BH, n = 1, with_ties = FALSE)

text_summary <- tibble(
  Item = c("AUC (ANOVA)", "AUC vs commercial (Tukey)", "AUC (sensitivity)",
           "Peak", "Days positive", "Daily prevalence (Fisher)"),
  Sentence = c(
    paste0("One-way ANOVA, F", auc_F$Df, ",", fit_auc$df.residual,
           " = ", sprintf("%.2f", auc_F$`F value`), ", ", fmt_p_text(auc_F$`Pr(>F)`), "."),
    paste0(vs_comm(auc_tukey, "p.value"), "."),
    paste0(kw_txt("auc_log10"), "; ", dunn_txt("auc_log10"), "."),
    paste0(kw_txt("peak_log10"), "; ", dunn_txt("peak_log10"), "."),
    paste0(kw_txt("viremic_days"), "; ", dunn_txt("viremic_days"), "."),
    paste0("Smallest BH-adjusted p = ", sprintf("%.3f", min_fisher$p_BH),
           " (day ", min_fisher$Day, ", ", min_fisher$group1, " ", min_fisher$pos1,
           " vs ", min_fisher$group2, " ", min_fisher$pos2, ").")
  )
)

cat("\n================ RESUMEN PARA EL TEXTO =================\n")
for (i in seq_len(nrow(text_summary))) {
  cat("\n", text_summary$Item[i], ": ", text_summary$Sentence[i], sep = "")
}
cat("\n\nTabla 3 (OR, IC 95 % y p con ajuste Tukey):\n")
print(or_table, n = Inf)


# ============================================================
# 19. FIGURA 5
#   A: cinética (media ± SEM), B: prevalencia (%),
#   C: AUC, D: pico, E: días positivos
# ============================================================

x_breaks <- time_levels

fig5A <- ggplot() +
  geom_point(
    data = datos,
    aes(x = Tiempo, y = log_load, color = Grupo, fill = Grupo, shape = Grupo),
    position = position_jitter(width = 0.25, height = 0, seed = 1),
    size = 1.2, alpha = 0.45
  ) +
  geom_line(data = viremia_summary,
            aes(x = Tiempo, y = mean_log, color = Grupo, group = Grupo),
            linewidth = 0.8) +
  geom_errorbar(data = viremia_summary,
                aes(x = Tiempo, ymin = pmax(0, mean_log - sem_log),
                    ymax = mean_log + sem_log, color = Grupo),
                width = 0.4, linewidth = 0.4) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = 0:7) +
  labs(x = "Days after challenge",
       y = expression("Viral load (" * log[10] * " copies/mL + 1)")) +
  theme_nature(10)

fig5B <- ggplot(viremia_summary,
                aes(x = Tiempo, y = pct_viremic, color = Grupo, group = Grupo)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Grupo, fill = Grupo), size = 2.4) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 105),
                     expand = c(0, 0), labels = function(x) paste0(x, "%")) +
  labs(x = "Days after challenge", y = "Viremic animals (%)") +
  theme_nature(10)

make_metric_boxplot <- function(metric, y_title, stat_df) {
  
  df_plot <- animal_metrics %>% select(Grupo, Value = all_of(metric))
  
  stat_sig <- stat_df %>% filter(p.adj < alpha_level)
  
  p <- ggplot(df_plot, aes(x = Grupo, y = Value)) +
    geom_boxplot(aes(fill = Grupo, color = Grupo), width = 0.45,
                 outlier.shape = NA, alpha = 0.15, linewidth = 0.4) +
    geom_point(aes(color = Grupo, fill = Grupo, shape = Grupo),
               position = position_jitter(width = 0.1, height = 0, seed = 1),
               size = 1.5, alpha = 0.85) +
    scale_x_discrete(labels = group_labels) +
    scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
    scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
    scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
    labs(x = NULL, y = y_title) +
    theme_nature(9) +
    theme(legend.position = "none")
  
  if (nrow(stat_sig) > 0) {
    stat_sig <- stat_sig %>%
      rstatix::add_y_position(data = animal_metrics,
                              formula = as.formula(paste(metric, "~ Grupo")),
                              step.increase = 0.12) %>%
      mutate(label = fmt_p_text(p.adj))
    
    p <- p + ggpubr::stat_pvalue_manual(stat_sig, label = "label",
                                        tip.length = 0.01, bracket.size = 0.35,
                                        size = 2.4)
  }
  
  p
}

fig5C <- make_metric_boxplot(
  "auc_log10",
  expression("AUC (" * log[10] * " copies/mL × days)"),
  auc_tukey_plot
)

fig5D <- make_metric_boxplot(
  "peak_log10",
  expression("Peak load (" * log[10] * " copies/mL)"),
  dunn_for_plot("peak_log10")
)

fig5E <- make_metric_boxplot(
  "viremic_days",
  "Days positive (n)",
  dunn_for_plot("viremic_days")
)

fig5_layout <- c(
  patchwork::area(1, 1, 3, 2),
  patchwork::area(4, 1, 6, 2),
  patchwork::area(1, 3, 2, 3),
  patchwork::area(3, 3, 4, 3),
  patchwork::area(5, 3, 6, 3)
)

fig5_final <- fig5A + fig5B + fig5C + fig5D + fig5E +
  plot_layout(design = fig5_layout, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top",
        plot.tag = element_text(face = "bold", size = 11))

save_nature(fig5_final, "Viremia_Figure_5", width = 9.0, height = 7.4)


# ============================================================
# 20. GUARDAR MODELOS
# ============================================================

saveRDS(list(auc_anova = fit_auc, prevalence_bglmer = fit_glmm),
        file.path(models_dir, "Viremia_modelos.rds"))


# ============================================================
# 21. INFORMACIÓN DEL ANÁLISIS
# ============================================================

analysis_info <- tibble(
  Item = c("Analysis", "Input file", "Number of animals", "Groups",
           "Days after challenge", "Positive sample", "Viral load scale",
           "AUC", "Peak and days positive", "Prevalence model",
           "Prevalence contrasts", "Daily prevalence", "R version"),
  Value = c(
    "Viremia after challenge (serum qPCR)",
    basename(input_file),
    as.character(n_distinct(datos$Animal)),
    paste(group_levels, collapse = ", "),
    paste(time_levels, collapse = ", "),
    if (is.null(positivity_column)) "Copias_mL > 0 (samples < LOQ counted as negative)"
    else paste0("Column '", positivity_column, "' (amplification detected)"),
    "log10(copies/mL + 1); samples < LOQ = 0 copies/mL",
    "Trapezoidal, days 0–21; one-way ANOVA + Tukey (primary); Kruskal–Wallis + Dunn-BH (sensitivity)",
    "Kruskal–Wallis + Dunn's test (BH)",
    paste0("blme::bglmer, binomial; group * day (categorical) + (1 | animal); days ",
           paste(glmm_days, collapse = ", "), "; normal prior SD ", prior_sd,
           " on fixed effects; gamma prior on random-effect covariance"),
    "Pairwise odds ratios averaged over days; Tukey-adjusted 95% CI and p (emmeans)",
    paste0("Pairwise Fisher's exact tests, BH within day; days ",
           paste(fisher_days, collapse = ", ")),
    R.version.string
  )
)


# ============================================================
# 22. EXPORTAR RESULTADOS A UN SOLO EXCEL
# ============================================================

output_excel <- file.path(tables_dir, "Viremia_resultados.xlsx")

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
  Sheet = c("Datos", "Resumen_dia", "Metricas_animal", "Metricas_resumen",
            "AUC_ANOVA", "AUC_supuestos", "AUC_Tukey", "Kruskal_Wallis",
            "Dunn_BH", "Tabla3_OR", "Fisher_dia", "Texto", "Analysis_info"),
  Description = c(
    "Processed data",
    "Viral load (log10) and prevalence by group and day",
    "AUC, peak, days positive and time to peak per animal",
    "Mean ± SD of the metrics by group",
    "One-way ANOVA of AUC",
    "Shapiro–Wilk and Levene tests for the AUC model",
    "Tukey pairwise comparisons of AUC",
    "Kruskal–Wallis tests (AUC, peak, days positive)",
    "Dunn's pairwise tests with BH adjustment",
    "Odds ratios from the penalized binomial GLMM (Table 3)",
    "Daily prevalence: pairwise Fisher's exact tests with BH adjustment",
    "Draft sentences for Section 3.4",
    "Analysis settings"
  )
)

add_sheet(wb, "README",           readme_table)
add_sheet(wb, "Datos",            datos %>% select(Animal, Grupo, Tiempo, Copias_mL, log_load, Viremic))
add_sheet(wb, "Resumen_dia",      viremia_summary)
add_sheet(wb, "Metricas_animal",  animal_metrics)
add_sheet(wb, "Metricas_resumen", metric_summary)
add_sheet(wb, "AUC_ANOVA",        auc_anova)
add_sheet(wb, "AUC_supuestos",    auc_assumptions)
add_sheet(wb, "AUC_Tukey",        auc_tukey)
add_sheet(wb, "Kruskal_Wallis",   kw_table)
add_sheet(wb, "Dunn_BH",          dunn_table)
add_sheet(wb, "Tabla3_OR",        or_table)
add_sheet(wb, "Fisher_dia",       fisher_by_day)
add_sheet(wb, "Texto",            text_summary)
add_sheet(wb, "Analysis_info",    analysis_info)

saveWorkbook(wb, output_excel, overwrite = TRUE)


# ============================================================
# 23. INFORME EN UN SOLO WORD
# ============================================================

output_word <- file.path(reports_dir, "Viremia_informe.docx")

make_ft <- function(df) {
  flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    autofit()
}

doc <- read_docx() %>%
  body_add_par("Viremia after challenge", style = "heading 1") %>%
  body_add_par(paste0("Generated on ", format(Sys.Date(), "%Y-%m-%d"),
                      " with ", R.version.string, "."), style = "Normal") %>%
  body_add_par("Analysis", style = "heading 2") %>%
  body_add_par(paste(
    "Viral load was analyzed as log10(copies/mL + 1). Cumulative exposure",
    "(trapezoidal AUC, days 0–21) was compared with one-way ANOVA and Tukey's",
    "test, with the Kruskal–Wallis test as a sensitivity analysis. Peak viral",
    "load and the number of positive sampling days were compared with the",
    "Kruskal–Wallis test followed by Dunn's test with Benjamini–Hochberg",
    "adjustment. The prevalence of viremia (days 3, 6 and 14) was analyzed",
    "with a penalized binomial generalized linear mixed model (blme::bglmer)",
    "including group, day and their interaction as fixed effects and animal as",
    "a random intercept; pairwise odds ratios averaged over days are reported",
    "with Tukey-adjusted 95% confidence intervals."
  ), style = "Normal") %>%
  body_add_par("Draft text for Section 3.4", style = "heading 2")

for (i in seq_len(nrow(text_summary))) {
  doc <- body_add_par(doc, paste0(text_summary$Item[i], ": ", text_summary$Sentence[i]),
                      style = "Normal")
}

doc <- doc %>%
  body_add_par("Metrics by group (mean ± SD)", style = "heading 2") %>%
  body_add_flextable(make_ft(
    metric_summary %>%
      transmute(Group = group_labels[as.character(Grupo)], n,
                `AUC (log10 copies/mL × days)` = sprintf("%.2f ± %.2f", auc_mean, auc_sd),
                `Peak (log10 copies/mL)`       = sprintf("%.2f ± %.2f", peak_log10_mean, peak_log10_sd),
                `Peak (copies/mL)`             = sprintf("%s ± %s",
                                                         format(round(peak_copies_mean), big.mark = ","),
                                                         format(round(peak_copies_sd), big.mark = ",")),
                `Days positive`                = sprintf("%.2f ± %.2f", viremic_days_mean, viremic_days_sd))
  )) %>%
  body_add_par("AUC: ANOVA assumptions", style = "heading 3") %>%
  body_add_flextable(make_ft(auc_assumptions %>% mutate(p_value = fmt_p_table(p_value)))) %>%
  body_add_par("Kruskal–Wallis tests", style = "heading 3") %>%
  body_add_flextable(make_ft(kw_table %>% select(Metric, n, H, df, p = p_text))) %>%
  body_add_par("Table 3. Pairwise odds ratios of viremia (days 3, 6 and 14)", style = "heading 2") %>%
  body_add_flextable(make_ft(or_table %>% select(Comparison, OR, `95% CI` = CI_95, `Adjusted p` = p_text))) %>%
  body_add_break() %>%
  body_add_par("Figure 5", style = "heading 2") %>%
  body_add_img(src = file.path(figures_dir, "Viremia_Figure_5.png"),
               width = 6.5, height = 6.5 * 7.4 / 9.0) %>%
  body_add_par(paste(
    "Viremia kinetics after challenge. (A) Serum viral load (log10 copies/mL + 1);",
    "lines and error bars show group means and SEM. (B) Prevalence of viremia.",
    "(C) Cumulative exposure (AUC), (D) peak viral load and (E) number of",
    "sampling days with detectable viremia. Boxplots show median and IQR.",
    "Brackets show adjusted p < 0.05 (Tukey for AUC; Dunn-BH for peak and days)."
  ), style = "Normal")

print(doc, target = output_word)


# ============================================================
# 24. SESSION INFO (reproducibilidad en GitHub)
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_viremia.txt"))


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
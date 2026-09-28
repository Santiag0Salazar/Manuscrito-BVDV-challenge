# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# TEMPERATURA RECTAL Y LINFOCITOS
#
# Secciones 3.1 (tolerancia posvacunal) y 3.6 (desafío),
# Figura 8A–B, Tabla 4 (filas de temperatura), Figuras S1 y S3
#
# Hojas del Excel:
#   "T° rectal desafio" : días 0–27 post-desafío (am y pm en días 0–2)
#   "T°Rectal v1"       : días 0–6 tras la primera dosis
#   "T°Rectal v2"       : días 0–2 tras el refuerzo
#   "LYM"               : linfocitos, promedios por grupo (sin datos individuales)
#
# Análisis (Secciones 2.10 y 2.11 del manuscrito):
#   - Temperatura: modelo lineal mixto (lmerTest::lmer), grupo * tiempo
#     (categórico), intercepto aleatorio por animal; F tipo III, Satterthwaite
#   - Sensibilidad: mismo modelo con correlación temporal CAR(1) (nlme)
#   - Día 8 (pico en control y comercial): ANOVA de una vía + Tukey
#   - Incidencia por animal (Tabla 4): ≥ 39.2 °C, > 40 °C, > 40 °C el día 8,
#     fiebre sostenida (≥ 2 lecturas consecutivas > 40 °C); Fisher exacto
#     recombinantes (agrupados) vs comercial
#   - Posvacunación: animales con > 40 °C y con fiebre sostenida por dosis
#   - Linfocitos: descriptivo (solo hay promedios por grupo)
#
# Salidas (carpeta results/ del proyecto):
#   figures/Temperature_Figure_8A_kinetics.(pdf|tiff|png)
#   figures/Temperature_Figure_8B_day8.(pdf|tiff|png)
#   figures/Temperature_Figure_S1_post_vaccination.(pdf|tiff|png)
#   figures/Lymphocytes_Figure_S3.(pdf|tiff|png)
#   tables/Temperatura_linfocitos_resultados.xlsx
#   reports/Temperatura_linfocitos_informe.docx
#   models/Temperatura_modelos.rds
#   sessionInfo_temperatura.txt
#
# Autor: Santiago Salazar
# ============================================================


# ============================================================
# 1. OPCIONES Y PAQUETES
# ============================================================

options(stringsAsFactors = FALSE)

required_packages <- c(
  "ragg",
  "systemfonts",
  "ggplot2", "dplyr", "tidyr", "tibble", "readxl", "openxlsx",
  "lme4", "lmerTest", "nlme", "emmeans", "patchwork", "cowplot",
  "purrr", "scales", "officer", "flextable"
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
  library(patchwork)
  library(cowplot)
  library(purrr)
  library(scales)
  library(officer)
  library(flextable)
})

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

# Tablas limpias generadas por 00_Limpieza_datos_clinicos.R
input_file     <- file.path(data_dir, "Temperaturas.xlsx")
input_file_lym <- file.path(data_dir, "Linfocitos.xlsx")

fever_threshold      <- 39.2   # Sección 2.10
high_fever_threshold <- 40.0   # categoría 3 de la Tabla 2 (> 40.0 °C)
day_peak             <- 8      # pico en los grupos control y comercial

# Los animales retirados y el cambio de rótulo 5451 -> 7766 (primera dosis)
# se aplican en 00_Limpieza_datos_clinicos.R
# (se repiten aquí solo para el registro del informe)
withdrawn_animals <- c("5451", "5452", "5434", "4906", "5442")
relabel_v1        <- c("5451" = "7766")

group_levels <- c("Control", "25ug", "50ug", "100ug", "Commercial")
recombinant  <- c("25ug", "50ug", "100ug")


# ============================================================
# 6. LECTURA Y LIMPIEZA
# ============================================================

for (f in c(input_file, input_file_lym)) {
  if (!file.exists(f)) stop(paste0("No se encontró el archivo:\n", f,
                                   "\nCorre primero 00_Limpieza_datos_clinicos.R"))
}

temperatures <- readxl::read_excel(input_file, sheet = "Temperaturas") %>%
  mutate(
    Animal = as.character(Animal),
    Grupo  = factor(Grupo, levels = group_levels),
    Day    = as.numeric(Day),
    Time   = as.numeric(Time),
    Temp   = as.numeric(Temp)
  )

if (any(is.na(temperatures$Grupo))) stop("Grupos no reconocidos en ", basename(input_file))

get_phase <- function(phase) {
  temperatures %>%
    filter(Phase == phase) %>%
    select(Animal, Grupo, Reading, Temp, Day, Session, Time) %>%
    mutate(Animal = factor(Animal)) %>%
    arrange(Grupo, Animal, Time)
}

temp_challenge <- get_phase("Challenge")
temp_v1        <- get_phase("First dose")
temp_v2        <- get_phase("Booster")

# Linfocitos: promedios por grupo
lym_long <- readxl::read_excel(input_file_lym, sheet = "Linfocitos") %>%
  arrange(Sampling_order) %>%
  mutate(
    Grupo = factor(Grupo, levels = group_levels),
    LYM   = as.numeric(LYM),
    Day   = as.numeric(Day)
  )

lym_head <- unique(lym_long$Sampling)
lym_long <- lym_long %>%
  mutate(Sampling = factor(Sampling, levels = lym_head)) %>%
  select(Grupo, Sampling, LYM, Phase, Day)


# ============================================================
# 7. COMPROBACIÓN BÁSICA
# ============================================================

cat("\n================ DATOS =================\n")
cat("\nAnimales por grupo (desafío):\n")
print(temp_challenge %>% distinct(Animal, Grupo) %>% count(Grupo))
cat("\nAnimales por grupo (primera dosis):\n")
print(temp_v1 %>% distinct(Animal, Grupo) %>% count(Grupo))
cat("\nAnimales por grupo (refuerzo):\n")
print(temp_v2 %>% distinct(Animal, Grupo) %>% count(Grupo))


# ============================================================
# 8. PALETA, FORMAS Y TEMA (iguales al resto de los scripts)
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

# Tipo de letra de todas las figuras (texto de ejes, leyendas y etiquetas)
base_family <- "Arial"

if (!base_family %in% systemfonts::system_fonts()$family) {
  warning("No se encontró la fuente ", base_family, "; se usará 'sans'.")
  base_family <- "sans"
}

# geom_text/annotate no heredan la fuente del tema
update_geom_defaults("text",  list(family = base_family))
update_geom_defaults("label", list(family = base_family))

theme_nature <- function(base_size = 10) {

  cowplot::theme_cowplot(font_size = base_size, font_family = base_family) +
    theme(
      plot.title           = element_text(face = "plain", size = base_size, hjust = 0.5),
      axis.title           = element_text(size = base_size),
      axis.text            = element_text(size = base_size, color = "grey15"),
      strip.background     = element_rect(fill = "grey95", color = NA),
      strip.text           = element_text(face = "bold", size = base_size),
      legend.title         = element_blank(),
      legend.position      = "top",
      legend.justification = "center",
      panel.grid.major.y   = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank()
    )
}


# ============================================================
# 9. FUNCIÓN PARA GUARDAR FIGURAS
# ============================================================

save_nature <- function(plot, filename, width = 7, height = 4) {

  # PDF con cairo_pdf y TIFF/PNG con ragg: usan la fuente Arial instalada en el sistema
  # y admiten caracteres Unicode (µ, °, ×)
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
# 10. FORMATO DE VALORES P
# ============================================================

fmt_p_table <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))
fmt_p_text  <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))


# ============================================================
# 11. RESUMEN POR GRUPO Y LECTURA
# ============================================================

summarise_temp <- function(df) {
  df %>%
    group_by(Grupo, Time, Day, Session) %>%
    summarise(
      n    = n(),
      mean = mean(Temp),
      sd   = sd(Temp),
      sem  = sd / sqrt(n),
      max  = max(Temp),
      .groups = "drop"
    )
}

summary_challenge <- summarise_temp(temp_challenge)
summary_v1        <- summarise_temp(temp_v1)
summary_v2        <- summarise_temp(temp_v2)


# ============================================================
# 12. MODELO LINEAL MIXTO: TEMPERATURA DURANTE EL DESAFÍO
# ============================================================

fit_lmm <- function(df) {
  df <- df %>% mutate(TimeF = factor(Time))
  lmerTest::lmer(Temp ~ Grupo * TimeF + (1 | Animal), data = df, REML = TRUE)
}

anova_lmm <- function(m, label) {
  a <- as.data.frame(anova(m, type = 3, ddf = "Satterthwaite"))
  tibble(
    Dataset = label,
    Effect  = recode(rownames(a), "Grupo" = "Group", "TimeF" = "Time",
                     "Grupo:TimeF" = "Group × time"),
    NumDF   = a$NumDF,
    DenDF   = round(a$DenDF, 1),
    F_value = round(a$`F value`, 2),
    p_value = a$`Pr(>F)`,
    p_text  = fmt_p_table(a$`Pr(>F)`)
  )
}

lmm_challenge <- fit_lmm(temp_challenge)
lmm_v1        <- fit_lmm(temp_v1)
lmm_v2        <- fit_lmm(temp_v2)

anova_table <- bind_rows(
  anova_lmm(lmm_challenge, "Challenge"),
  anova_lmm(lmm_v1,        "First dose"),
  anova_lmm(lmm_v2,        "Booster")
)

# Sensibilidad: correlación temporal CAR(1) entre lecturas del mismo animal
lmm_car1 <- tryCatch(
  nlme::lme(
    Temp ~ Grupo * TimeF,
    random      = ~ 1 | Animal,
    correlation = nlme::corCAR1(form = ~ Time | Animal),
    data        = temp_challenge %>% mutate(TimeF = factor(Time)),
    method      = "REML",
    control     = nlme::lmeControl(opt = "optim")
  ),
  error = function(e) { message("CAR(1) no convergió: ", e$message); NULL }
)

anova_car1 <- if (!is.null(lmm_car1)) {
  a <- as.data.frame(anova(lmm_car1, type = "marginal"))
  tibble(
    Effect  = rownames(a), numDF = a$numDF, denDF = a$denDF,
    F_value = round(a$`F-value`, 2), p_value = a$`p-value`,
    p_text  = fmt_p_table(a$`p-value`),
    Phi     = round(as.numeric(coef(lmm_car1$modelStruct$corStruct,
                                    unconstrained = FALSE)), 3)
  ) %>% filter(Effect != "(Intercept)")
} else tibble()


# ============================================================
# 13. DÍA 8: ANOVA DE UNA VÍA + TUKEY
# ============================================================

temp_day8 <- temp_challenge %>% filter(Time == day_peak)

fit_day8 <- lm(Temp ~ Grupo, data = temp_day8)

day8_anova <- as.data.frame(anova(fit_day8)) %>%
  rownames_to_column("Effect") %>%
  as_tibble()

day8_tukey <- as.data.frame(pairs(emmeans(fit_day8, ~ Grupo), adjust = "tukey")) %>%
  as_tibble() %>%
  mutate(p_text = fmt_p_table(p.value))

day8_descriptives <- temp_day8 %>%
  group_by(Grupo) %>%
  summarise(n = n(), mean = mean(Temp), sd = sd(Temp),
            min = min(Temp), max = max(Temp), .groups = "drop")


# ============================================================
# 14. INCIDENCIA POR ANIMAL (TABLA 4) Y FIEBRE SOSTENIDA
# ============================================================

# Fiebre sostenida: dos o más lecturas consecutivas sobre el umbral
max_run <- function(x) {
  r <- rle(x)
  if (any(r$values)) max(r$lengths[r$values]) else 0L
}

animal_outcomes <- function(df) {
  df %>%
    arrange(Animal, Time) %>%
    group_by(Grupo, Animal) %>%
    summarise(
      max_temp            = max(Temp),
      any_ge_fever_am     = any(Temp[Session == "am"] >= fever_threshold),
      any_ge_fever_all    = any(Temp >= fever_threshold),
      n_am_ge_fever       = sum(Temp[Session == "am"] >= fever_threshold),
      any_gt_40           = any(Temp > high_fever_threshold),
      gt_40_day_peak      = any(Temp[Time == day_peak] > high_fever_threshold),
      ge_fever_day_peak   = any(Temp[Time == day_peak] >= fever_threshold),
      sustained_gt_40     = max_run(Temp > high_fever_threshold) >= 2,
      .groups = "drop"
    )
}

outcomes_challenge <- animal_outcomes(temp_challenge)
outcomes_v1        <- animal_outcomes(temp_v1)
outcomes_v2        <- animal_outcomes(temp_v2)

count_by_group <- function(outcomes, label) {
  outcomes %>%
    group_by(Grupo) %>%
    summarise(
      n = n(),
      across(c(any_ge_fever_am, any_ge_fever_all, any_gt_40,
               gt_40_day_peak, ge_fever_day_peak, sustained_gt_40), sum),
      .groups = "drop"
    ) %>%
    mutate(Dataset = label, .before = 1)
}

incidence_table <- bind_rows(
  count_by_group(outcomes_challenge, "Challenge"),
  count_by_group(outcomes_v1,        "First dose"),
  count_by_group(outcomes_v2,        "Booster")
)

# Fisher: recombinantes (agrupados) vs comercial
fisher_rec_vs_comm <- function(outcomes, variable, label) {

  sub <- outcomes %>%
    filter(Grupo %in% c(recombinant, "Commercial")) %>%
    mutate(Set = ifelse(Grupo == "Commercial", "Commercial", "Recombinant"))

  tab <- table(factor(sub$Set, levels = c("Recombinant", "Commercial")),
               factor(sub[[variable]], levels = c(FALSE, TRUE)))

  tibble(
    Dataset     = label,
    Outcome     = variable,
    Recombinant = paste0(tab["Recombinant", "TRUE"], "/", sum(tab["Recombinant", ])),
    Commercial  = paste0(tab["Commercial", "TRUE"], "/", sum(tab["Commercial", ])),
    p_value     = fisher.test(tab)$p.value,
    p_text      = fmt_p_table(fisher.test(tab)$p.value)
  )
}

fisher_table <- bind_rows(
  fisher_rec_vs_comm(outcomes_challenge, "ge_fever_day_peak", "Challenge"),
  fisher_rec_vs_comm(outcomes_challenge, "any_gt_40",         "Challenge"),
  fisher_rec_vs_comm(outcomes_challenge, "any_ge_fever_am",   "Challenge"),
  fisher_rec_vs_comm(outcomes_challenge, "sustained_gt_40",   "Challenge"),
  fisher_rec_vs_comm(outcomes_v1,        "any_gt_40",         "First dose"),
  fisher_rec_vs_comm(outcomes_v2,        "any_gt_40",         "Booster"),
  fisher_rec_vs_comm(outcomes_v2,        "sustained_gt_40",   "Booster")
)


# ============================================================
# 15. LINFOCITOS: CAMBIO RESPECTO DEL DÍA DEL DESAFÍO
# ============================================================

lym_challenge <- lym_long %>%
  filter(Phase == "Challenge") %>%
  group_by(Grupo) %>%
  mutate(
    baseline   = LYM[Day == 0][1],
    pct_change = 100 * (LYM - baseline) / baseline
  ) %>%
  ungroup()


# ============================================================
# 16. FRASES LISTAS PARA EL MANUSCRITO
# ============================================================

fmt_F <- function(df, effect) {
  r <- df %>% filter(Effect == effect)
  paste0("F", r$NumDF, ",", round(r$DenDF), " = ", sprintf("%.2f", r$F_value),
         ", ", fmt_p_text(r$p_value))
}

a_ch  <- anova_table %>% filter(Dataset == "Challenge")
d8_F  <- day8_anova %>% filter(Effect == "Grupo")
ctrl8 <- day8_descriptives %>% filter(Grupo == "Control")

tukey_vs_comm <- day8_tukey %>%
  filter(grepl("Commercial", contrast), !grepl("Control", contrast)) %>%
  mutate(txt = paste0(gsub("ug", " µg", contrast), ", ", fmt_p_text(p.value))) %>%
  pull(txt) %>% paste(collapse = "; ")

tukey_among_rec <- day8_tukey %>%
  filter(!grepl("Commercial|Control", contrast))

f_row <- function(label, outcome) fisher_table %>% filter(Dataset == label, Outcome == outcome)

text_summary <- tibble(
  Item = c("Temperature LMM", "Day 8", "Day 8 control",
           "Fever ≥ 39.2 °C on day 8", "Temperature > 40 °C, whole period",
           "Post-vaccination, first dose", "Post-vaccination, booster"),
  Sentence = c(
    paste0("Time, ", fmt_F(a_ch, "Time"), "; group × time, ", fmt_F(a_ch, "Group × time"),
           "; group, ", fmt_F(a_ch, "Group"), "."),
    paste0("One-way ANOVA, F", d8_F$Df, ",", fit_day8$df.residual, " = ",
           sprintf("%.2f", d8_F$`F value`), ", ", fmt_p_text(d8_F$`Pr(>F)`),
           ". Tukey vs commercial: ", tukey_vs_comm,
           ". Among recombinant doses: all ", ifelse(all(tukey_among_rec$p.value > 0.9),
                                                    "p > 0.9", paste0("p ≥ ", sprintf("%.3f", min(tukey_among_rec$p.value)))),
           "."),
    sprintf("Control (n = %d): %.2f ± %.2f °C (mean ± SD).", ctrl8$n, ctrl8$mean, ctrl8$sd),
    with(f_row("Challenge", "ge_fever_day_peak"),
         paste0("Recombinant ", Recombinant, " vs commercial ", Commercial, " (Fisher, ", fmt_p_text(p_value), ").")),
    with(f_row("Challenge", "any_gt_40"),
         paste0("Recombinant ", Recombinant, " vs commercial ", Commercial, " (Fisher, ", fmt_p_text(p_value), ").")),
    with(f_row("First dose", "any_gt_40"),
         paste0("> 40 °C: recombinant ", Recombinant, " vs commercial ", Commercial, " (Fisher, ", fmt_p_text(p_value), ").")),
    paste0(
      with(f_row("Booster", "any_gt_40"),
           paste0("> 40 °C: recombinant ", Recombinant, " vs commercial ", Commercial, "; ")),
      with(f_row("Booster", "sustained_gt_40"),
           paste0("sustained fever: recombinant ", Recombinant, " vs commercial ", Commercial, ".")))
  )
)

cat("\n================ RESUMEN PARA EL TEXTO =================\n")
for (i in seq_len(nrow(text_summary))) {
  cat("\n", text_summary$Item[i], ": ", text_summary$Sentence[i], sep = "")
}
cat("\n\nIncidencia por grupo:\n")
print(incidence_table, n = Inf, width = Inf)


# ============================================================
# 17. FIGURA 8A: TEMPERATURA DURANTE EL DESAFÍO
# ============================================================

fig8A <- ggplot() +
  geom_line(data = temp_challenge,
            aes(x = Time, y = Temp, group = Animal, color = Grupo),
            linewidth = 0.3, alpha = 0.18) +
  geom_line(data = summary_challenge,
            aes(x = Time, y = mean, color = Grupo, group = Grupo),
            linewidth = 0.8) +
  geom_errorbar(data = summary_challenge,
                aes(x = Time, ymin = mean - sem, ymax = mean + sem, color = Grupo),
                width = 0.15, linewidth = 0.3) +
  geom_point(data = summary_challenge,
             aes(x = Time, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 1.8) +
  geom_hline(yintercept = fever_threshold, linetype = "dashed",
             linewidth = 0.45, color = "grey30") +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 28, 2), limits = c(-0.3, 28)) +
  scale_y_continuous(breaks = seq(37.5, 42, 0.5)) +
  labs(x = "Days after challenge", y = "Rectal temperature (°C)") +
  theme_nature(10)

save_nature(fig8A, "Temperature_Figure_8A_kinetics", width = 6.5, height = 3.6)


# ============================================================
# 18. FIGURA 8B: DÍA 8 (BARRAS SOLO PARA p < 0.05, SIN EL CONTROL)
# ============================================================

x_pos <- setNames(seq_along(group_levels), group_levels)

brackets_day8 <- day8_tukey %>%
  filter(p.value < 0.05, !grepl("Control", contrast)) %>%
  mutate(
    g1 = trimws(sub(" - .*$", "", contrast)),
    g2 = trimws(sub("^.* - ", "", contrast)),
    x  = pmin(x_pos[g1], x_pos[g2]),
    xend = pmax(x_pos[g1], x_pos[g2])
  ) %>%
  arrange(xend - x, x)

y_top   <- max(temp_day8$Temp)
y_range <- diff(range(temp_day8$Temp))

brackets_day8 <- brackets_day8 %>%
  mutate(
    y      = y_top + 0.10 * y_range * row_number(),
    y_text = y + 0.02 * y_range,
    label  = fmt_p_text(p.value)
  )

fig8B <- ggplot(temp_day8, aes(x = Grupo, y = Temp)) +
  geom_boxplot(aes(color = Grupo, fill = Grupo), width = 0.5,
               outlier.shape = NA, alpha = 0.15, linewidth = 0.4) +
  geom_point(aes(color = Grupo, fill = Grupo, shape = Grupo),
             position = position_jitter(width = 0.1, height = 0, seed = 1),
             size = 1.6, alpha = 0.9) +
  geom_segment(data = brackets_day8, aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.35, color = "grey30") +
  geom_text(data = brackets_day8, aes(x = (x + xend) / 2, y = y_text, label = label),
            inherit.aes = FALSE, size = 2.3, vjust = 0) +
  geom_hline(yintercept = fever_threshold, linetype = "dashed",
             linewidth = 0.45, color = "grey30") +
  scale_x_discrete(labels = group_labels) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  labs(x = NULL, y = "Rectal temperature at day 8 (°C)") +
  guides(color = "none", fill = "none", shape = "none") +
  theme_nature(10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

save_nature(fig8B, "Temperature_Figure_8B_day8", width = 2.6, height = 3.8)

# Paneles 8A y 8B guardados para ensamblar la Figura 8 completa en Diarrea_post_desafio.R
saveRDS(list(fig8A = fig8A, fig8B = fig8B),
        file.path(models_dir, "Temperature_Figure_8AB_plots.rds"))


# ============================================================
# 19. FIGURA S1: TEMPERATURA POSVACUNAL (AMBAS DOSIS)
# ============================================================

post_vac <- bind_rows(
  summary_v1 %>% mutate(Dose = "First dose"),
  summary_v2 %>% mutate(Dose = "Booster")
) %>%
  mutate(Dose = factor(Dose, levels = c("First dose", "Booster")))

post_vac_raw <- bind_rows(
  temp_v1 %>% mutate(Dose = "First dose"),
  temp_v2 %>% mutate(Dose = "Booster")
) %>%
  mutate(Dose = factor(Dose, levels = c("First dose", "Booster")))

figS1 <- ggplot() +
  geom_line(data = post_vac_raw,
            aes(x = Time, y = Temp, group = Animal, color = Grupo),
            linewidth = 0.3, alpha = 0.18) +
  geom_line(data = post_vac, aes(x = Time, y = mean, color = Grupo, group = Grupo),
            linewidth = 0.8) +
  geom_errorbar(data = post_vac,
                aes(x = Time, ymin = mean - sem, ymax = mean + sem, color = Grupo),
                width = 0.1, linewidth = 0.3) +
  geom_point(data = post_vac,
             aes(x = Time, y = mean, color = Grupo, fill = Grupo, shape = Grupo),
             size = 1.8) +
  geom_hline(yintercept = high_fever_threshold, linetype = "dashed",
             linewidth = 0.45, color = "grey30") +
  facet_wrap(~ Dose, scales = "free_x") +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Days after vaccination (evening readings at + 0.5)",
       y = "Rectal temperature (°C)") +
  theme_nature(10)

save_nature(figS1, "Temperature_Figure_S1_post_vaccination", width = 7.5, height = 3.6)


# ============================================================
# 20. FIGURA S3: LINFOCITOS (PROMEDIOS POR GRUPO)
# ============================================================

lym_plot <- lym_long %>%
  mutate(x = as.numeric(Sampling))

challenge_x <- lym_plot %>% filter(Phase == "Challenge", Day == 0) %>% pull(x) %>% unique()

x_labels <- lym_plot %>%
  distinct(x, Phase, Day) %>%
  arrange(x) %>%
  mutate(lab = ifelse(Phase == "Challenge",
                      paste0("C+", Day), paste0("V", Day)))

figS3 <- ggplot(lym_plot, aes(x = x, y = LYM, color = Grupo, group = Grupo)) +
  geom_vline(xintercept = challenge_x, linetype = "dashed",
             linewidth = 0.45, color = "grey30") +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Grupo, fill = Grupo), size = 2.2) +
  scale_x_continuous(breaks = x_labels$x, labels = x_labels$lab) +
  scale_color_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = palette_nature, labels = group_labels, drop = FALSE) +
  scale_shape_manual(values = shape_groups, labels = group_labels, drop = FALSE) +
  labs(x = "Sampling (V = days after first vaccination; C = days after challenge)",
       y = expression("Lymphocytes (" * 10^9 * "/L)")) +
  theme_nature(10)

save_nature(figS3, "Lymphocytes_Figure_S3", width = 7.0, height = 3.6)


# ============================================================
# 21. GUARDAR MODELOS
# ============================================================

saveRDS(list(challenge = lmm_challenge, first_dose = lmm_v1, booster = lmm_v2,
             challenge_car1 = lmm_car1, day8 = fit_day8),
        file.path(models_dir, "Temperatura_modelos.rds"))


# ============================================================
# 22. INFORMACIÓN DEL ANÁLISIS
# ============================================================

analysis_info <- tibble(
  Item = c("Input file", "Animals (challenge)", "Withdrawn animals excluded",
           "Relabelled rows", "Time coding", "Fever threshold",
           "High fever", "Sustained fever", "Temperature model",
           "Sensitivity model", "Day 8", "Incidence tests", "Lymphocytes",
           "R version"),
  Value = c(
    basename(input_file),
    as.character(n_distinct(temp_challenge$Animal)),
    paste(withdrawn_animals, collapse = ", "),
    paste0(names(relabel_v1), " -> ", relabel_v1, " (first-dose sheet)"),
    "Sheet 'Dia k' = day k − 1; evening readings = day + 0.5",
    paste0("≥ ", fever_threshold, " °C"),
    paste0("> ", high_fever_threshold, " °C"),
    paste0("≥ 2 consecutive readings > ", high_fever_threshold, " °C"),
    "lmerTest::lmer, group * time (categorical) + (1 | animal); type III F, Satterthwaite",
    "nlme::lme with CAR(1) within-animal correlation",
    paste0("Day ", day_peak, ": one-way ANOVA + Tukey (emmeans)"),
    "Fisher's exact test, recombinant groups pooled vs commercial",
    "Group means only (no individual values): descriptive",
    R.version.string
  )
)


# ============================================================
# 23. EXPORTAR RESULTADOS A UN SOLO EXCEL
# ============================================================

output_excel <- file.path(tables_dir, "Temperatura_linfocitos_resultados.xlsx")

wb <- createWorkbook()

add_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  if (is.null(data) || nrow(data) == 0) data <- data.frame(Information = "No data available")
  writeData(wb, sheet = sheet_name, x = data)
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(data)), widths = "auto")
}

add_sheet(wb, "Datos_desafio",      temp_challenge)
add_sheet(wb, "Datos_dosis1",       temp_v1)
add_sheet(wb, "Datos_refuerzo",     temp_v2)
add_sheet(wb, "Resumen_desafio",    summary_challenge)
add_sheet(wb, "Resumen_dosis1",     summary_v1)
add_sheet(wb, "Resumen_refuerzo",   summary_v2)
add_sheet(wb, "LMM_ANOVA",          anova_table)
add_sheet(wb, "LMM_CAR1",           anova_car1)
add_sheet(wb, "Dia8_ANOVA",         day8_anova)
add_sheet(wb, "Dia8_Tukey",         day8_tukey)
add_sheet(wb, "Dia8_descriptivos",  day8_descriptives)
add_sheet(wb, "Por_animal_desafio", outcomes_challenge)
add_sheet(wb, "Por_animal_dosis1",  outcomes_v1)
add_sheet(wb, "Por_animal_refuerzo", outcomes_v2)
add_sheet(wb, "Incidencia_Tabla4",  incidence_table)
add_sheet(wb, "Fisher",             fisher_table)
add_sheet(wb, "Linfocitos",         lym_long)
add_sheet(wb, "Linfocitos_desafio", lym_challenge)
add_sheet(wb, "Texto",              text_summary)
add_sheet(wb, "Analysis_info",      analysis_info)

saveWorkbook(wb, output_excel, overwrite = TRUE)


# ============================================================
# 24. INFORME EN UN SOLO WORD
# ============================================================

output_word <- file.path(reports_dir, "Temperatura_linfocitos_informe.docx")

make_ft <- function(df) {
  flextable(df) %>% theme_booktabs() %>% fontsize(size = 9, part = "all") %>% autofit()
}

doc <- read_docx() %>%
  body_add_par("Rectal temperature and lymphocytes", style = "heading 1") %>%
  body_add_par(paste0("Generated on ", format(Sys.Date(), "%Y-%m-%d"),
                      " with ", R.version.string, "."), style = "Normal") %>%
  body_add_par("Draft text", style = "heading 2")

for (i in seq_len(nrow(text_summary))) {
  doc <- body_add_par(doc, paste0(text_summary$Item[i], ": ", text_summary$Sentence[i]),
                      style = "Normal")
}

doc <- doc %>%
  body_add_par("Mixed models (temperature)", style = "heading 2") %>%
  body_add_flextable(make_ft(anova_table %>% select(Dataset, Effect, NumDF, DenDF, F_value, p = p_text))) %>%
  body_add_par("Sensitivity: CAR(1) within-animal correlation (challenge)", style = "heading 3") %>%
  body_add_flextable(make_ft(if (nrow(anova_car1) > 0)
    anova_car1 %>% select(Effect, numDF, denDF, F_value, p = p_text, Phi)
    else data.frame(Information = "Model did not converge"))) %>%
  body_add_par("Day 8: Tukey comparisons", style = "heading 2") %>%
  body_add_flextable(make_ft(day8_tukey %>%
                               transmute(Contrast = gsub("ug", " µg", contrast),
                                         Estimate = round(estimate, 2),
                                         t = round(t.ratio, 2), p = p_text))) %>%
  body_add_par("Incidence by group (Table 4 and post-vaccination)", style = "heading 2") %>%
  body_add_flextable(make_ft(incidence_table %>%
                               mutate(Grupo = group_labels[as.character(Grupo)]) %>%
                               rename(Group = Grupo,
                                      `≥39.2 am` = any_ge_fever_am,
                                      `≥39.2 all` = any_ge_fever_all,
                                      `>40` = any_gt_40,
                                      `>40 day 8` = gt_40_day_peak,
                                      `≥39.2 day 8` = ge_fever_day_peak,
                                      `Sustained >40` = sustained_gt_40))) %>%
  body_add_par("Fisher's exact tests (recombinant pooled vs commercial)", style = "heading 2") %>%
  body_add_flextable(make_ft(fisher_table %>% select(Dataset, Outcome, Recombinant, Commercial, p = p_text))) %>%
  body_add_break() %>%
  body_add_par("Figure 8A–B", style = "heading 2") %>%
  body_add_img(file.path(figures_dir, "Temperature_Figure_8A_kinetics.png"),
               width = 6.3, height = 6.3 * 3.6 / 6.5) %>%
  body_add_img(file.path(figures_dir, "Temperature_Figure_8B_day8.png"),
               width = 2.4, height = 2.4 * 3.8 / 2.6) %>%
  body_add_par("Figure S1: post-vaccination temperature", style = "heading 2") %>%
  body_add_img(file.path(figures_dir, "Temperature_Figure_S1_post_vaccination.png"),
               width = 6.3, height = 6.3 * 3.6 / 7.5) %>%
  body_add_par("Figure S3: lymphocytes (group means)", style = "heading 2") %>%
  body_add_img(file.path(figures_dir, "Lymphocytes_Figure_S3.png"),
               width = 6.3, height = 6.3 * 3.6 / 7.0)

print(doc, target = output_word)


# ============================================================
# 25. SESSION INFO Y MENSAJE FINAL
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo_temperatura.txt"))

cat("\n\n============================================================\n")
cat("ANÁLISIS COMPLETADO CORRECTAMENTE\n")
cat("============================================================\n")
cat("\nExcel:  ", output_excel, "\n")
cat("Word:   ", output_word, "\n")
cat("Figuras:", figures_dir, "\n")
cat("\n============================================================\n")

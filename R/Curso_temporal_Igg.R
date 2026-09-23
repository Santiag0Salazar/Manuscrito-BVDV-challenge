system("git add R/Curso_temporal_Igg.R")
system("git add R/Curso_temporal_Igg.R data/processed/Curso_temporal_Igg.xlsx")
system("git status")
system('git commit -m "Add processed IgG data and temporal analysis script"')
system("git push")

options(stringsAsFactors = FALSE)

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "nlme", "emmeans",
  "broom", "patchwork", "cowplot", "pracma", "scales", "splines"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(nlme)
  library(emmeans)
  library(broom)
  library(patchwork)
  library(cowplot)
  library(pracma)
  library(scales)
  library(splines)
})

set.seed(20260620)

base_dir <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
out_dir <- getwd()
tables_dir <- file.path(out_dir, "tablas")
figures_dir <- file.path(out_dir, "figuras")
models_dir <- file.path(out_dir, "modelos")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)

read_elisa_data <- function(path) {
  read.csv(path, sep = ";", dec = ".", check.names = FALSE) |>
    mutate(
      Animal = factor(Animal),
      Grupo = factor(Grupo, levels = c("Control", "25ug", "50ug", "100ug", "Commercial")),
      Tiempo = as.numeric(Tiempo),
      Absorbancia = as.numeric(Absorbancia),
      TiempoF = factor(Tiempo, levels = sort(unique(Tiempo))),
      TimeIndex = as.numeric(factor(Tiempo, levels = sort(unique(Tiempo))))
    ) |>
    arrange(Grupo, Animal, Tiempo)
}

datos <- read_elisa_data(file.path(base_dir, "Datos.csv"))
datos_vac <- datos |> filter(Grupo != "Control") |> droplevels()
tiempos <- sort(unique(datos$Tiempo))

# Paleta Nature para los puntos
palette_nature <- c(
  "Control"   = "#4D4D4D", # Gris oscuro
  "25ug"      = "#004B87", # Azul oscuro
  "50ug"      = "#E69F00", # Naranja
  "100ug"     = "#2ECC71", # Verde claro manzana
  "Comercial" = "#D55E00"  # Rojo bermellón
)

theme_nature <- function(base_size = 8) {
  cowplot::theme_cowplot(font_size = base_size, font_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size, color = "grey25"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, color = "grey15"),
      legend.title = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = base_size - 1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
}

save_nature <- function(plot, filename, width, height) {
  ggsave(file.path(figures_dir, paste0(filename, ".pdf")), plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figures_dir, paste0(filename, ".tiff")), plot, width = width, height = height, units = "in", dpi = 600, compression = "lzw")
  ggsave(file.path(figures_dir, paste0(filename, ".png")), plot, width = width, height = height, units = "in", dpi = 300)
}

mean_ci <- datos |>
  group_by(Grupo, Tiempo) |>
  summarise(
    n = n(),
    mean = mean(Absorbancia),
    sd = sd(Absorbancia),
    se = sd / sqrt(n),
    ci95 = qt(0.975, df = pmax(n - 1, 1)) * se,
    lower = mean - ci95,
    upper = mean + ci95,
    .groups = "drop"
  )

group_summary <- datos |>
  group_by(Grupo, Tiempo) |>
  summarise(
    n_animals = n_distinct(Animal),
    mean_abs = mean(Absorbancia),
    sd_abs = sd(Absorbancia),
    median_abs = median(Absorbancia),
    min_abs = min(Absorbancia),
    max_abs = max(Absorbancia),
    .groups = "drop"
  )
write_csv(group_summary, file.path(tables_dir, "01_resumen_por_grupo_y_tiempo.csv"))

animal_metrics <- datos |>
  group_by(Grupo, Animal) |>
  arrange(Tiempo, .by_group = TRUE) |>
  summarise(
    baseline_d0 = Absorbancia[Tiempo == 0][1],
    day14 = Absorbancia[Tiempo == 14][1],
    day21 = Absorbancia[Tiempo == 21][1],
    day28 = Absorbancia[Tiempo == 28][1],
    day56 = Absorbancia[Tiempo == 56][1],
    primary_peak_0_21 = max(Absorbancia[Tiempo <= 21]),
    post_peak_28_56 = max(Absorbancia[Tiempo >= 28]),
    total_peak = max(Absorbancia),
    time_to_peak = Tiempo[which.max(Absorbancia)][1],
    booster_delta_21_28 = day28 - day21,
    late_delta_28_56 = day56 - day28,
    persistence_ratio_56_peak = day56 / total_peak,
    auc_0_21 = pracma::trapz(Tiempo[Tiempo <= 21], Absorbancia[Tiempo <= 21]),
    auc_21_56 = pracma::trapz(Tiempo[Tiempo >= 21], Absorbancia[Tiempo >= 21]),
    auc_total = pracma::trapz(Tiempo, Absorbancia),
    .groups = "drop"
  )
write_csv(animal_metrics, file.path(tables_dir, "02_metricas_cineticas_por_animal.csv"))

fit_lme <- function(expr, data, method = "REML", correlation = NULL) {
  tryCatch(
    nlme::lme(
      fixed = expr,
      random = ~ 1 | Animal,
      data = data,
      method = method,
      correlation = correlation,
      control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
    ),
    error = function(e) NULL
  )
}

model_candidates <- list(
  measured_time = fit_lme(Absorbancia ~ Grupo * TiempoF, datos, "ML"),
  measured_time_ar1 = fit_lme(Absorbancia ~ Grupo * TiempoF, datos, "ML", corAR1(form = ~ TimeIndex | Animal)),
  spline_df3 = fit_lme(Absorbancia ~ Grupo * ns(Tiempo, df = 3), datos, "ML"),
  spline_df4 = fit_lme(Absorbancia ~ Grupo * ns(Tiempo, df = 4), datos, "ML"),
  cubic = fit_lme(Absorbancia ~ Grupo * poly(Tiempo, 3), datos, "ML")
)

model_comparison <- bind_rows(lapply(names(model_candidates), function(nm) {
  fit <- model_candidates[[nm]]
  if (is.null(fit)) {
    data.frame(model = nm, converged = FALSE, AIC = NA_real_, BIC = NA_real_, logLik = NA_real_)
  } else {
    data.frame(model = nm, converged = TRUE, AIC = AIC(fit), BIC = BIC(fit), logLik = as.numeric(logLik(fit)))
  }
})) |>
  arrange(AIC)
write_csv(model_comparison, file.path(tables_dir, "03_comparacion_modelos_mixtos_AIC.csv"))

final_time_model <- fit_lme(Absorbancia ~ Grupo * TiempoF, datos, "REML")
saveRDS(final_time_model, file.path(models_dir, "modelo_mixto_tiempo_medido_REML.rds"))

anova_time <- as.data.frame(anova(final_time_model)) |>
  tibble::rownames_to_column("term")
write_csv(anova_time, file.path(tables_dir, "04_anova_modelo_mixto_tiempo_medido.csv"))

emm_day <- emmeans(final_time_model, ~ Grupo | TiempoF, data = datos)
emm_day_table <- as.data.frame(summary(emm_day, infer = TRUE)) |>
  rename(Tiempo = TiempoF)
write_csv(emm_day_table, file.path(tables_dir, "05_medias_marginales_por_dia.csv"))

pairwise_day <- as.data.frame(summary(pairs(emm_day, adjust = "tukey"))) |>
  rename(Tiempo = TiempoF)
write_csv(pairwise_day, file.path(tables_dir, "06_comparaciones_tukey_por_dia.csv"))

metric_tests <- function(metric, data = animal_metrics |> filter(Grupo != "Control") |> droplevels()) {
  f <- reformulate("Grupo", response = metric)
  fit <- lm(f, data = data)
  emm <- emmeans(fit, ~ Grupo, data = data)
  list(
    anova = broom::tidy(anova(fit)) |> mutate(metric = metric, .before = 1),
    contrasts = as.data.frame(summary(pairs(emm, adjust = "tukey"))) |> mutate(metric = metric, .before = 1),
    means = as.data.frame(summary(emm, infer = TRUE)) |> mutate(metric = metric, .before = 1)
  )
}

metrics_to_test <- c(
  "primary_peak_0_21", "post_peak_28_56", "booster_delta_21_28",
  "late_delta_28_56", "persistence_ratio_56_peak", "auc_0_21",
  "auc_21_56", "auc_total"
)
metric_results <- lapply(metrics_to_test, metric_tests)
write_csv(bind_rows(lapply(metric_results, `[[`, "anova")), file.path(tables_dir, "07_anova_metricas_por_animal.csv"))
write_csv(bind_rows(lapply(metric_results, `[[`, "contrasts")), file.path(tables_dir, "08_tukey_metricas_por_animal.csv"))
write_csv(bind_rows(lapply(metric_results, `[[`, "means")), file.path(tables_dir, "09_medias_marginales_metricas.csv"))

dose_trend_data <- animal_metrics |>
  filter(Grupo %in% c("25ug", "50ug", "100ug")) |>
  mutate(Dose = as.numeric(sub("ug", "", as.character(Grupo))))
dose_trends <- bind_rows(lapply(metrics_to_test, function(metric) {
  fit <- lm(reformulate("Dose", response = metric), data = dose_trend_data)
  broom::tidy(fit) |>
    filter(term == "Dose") |>
    mutate(metric = metric, .before = 1, r_squared = summary(fit)$r.squared)
}))
write_csv(dose_trends, file.path(tables_dir, "10_tendencia_dosis_recombinantes.csv"))

group_means <- datos_vac |>
  group_by(Grupo, Tiempo) |>
  summarise(Absorbancia = mean(Absorbancia), .groups = "drop")

fit_nls_models <- function(df) {
  y_min <- min(df$Absorbancia)
  y_max <- max(df$Absorbancia)
  amp <- max(y_max - y_min, 0.1)
  safe_nls <- function(formula, start, lower, upper, model_name) {
    fit <- tryCatch(
      suppressWarnings(nls(
        formula,
        data = df,
        start = start,
        algorithm = "port",
        lower = lower,
        upper = upper,
        control = nls.control(maxiter = 500, warnOnly = TRUE)
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    pred <- predict(fit, df)
    list(
      model = model_name,
      fit = fit,
      params = coef(fit),
      AIC = AIC(fit),
      RSS = sum(resid(fit)^2),
      RMSE = sqrt(mean(resid(fit)^2)),
      R2 = 1 - sum((df$Absorbancia - pred)^2) / sum((df$Absorbancia - mean(df$Absorbancia))^2)
    )
  }

  list(
    Logistic = safe_nls(
      Absorbancia ~ Bottom + (Top - Bottom) / (1 + exp(-(Tiempo - T50) / s)),
      start = list(Bottom = y_min, Top = y_max, T50 = 18, s = 3),
      lower = c(Bottom = 0, Top = 0.3, T50 = 0, s = 0.2),
      upper = c(Bottom = 1.0, Top = 3.0, T50 = 56, s = 20),
      model_name = "Logistic"
    ),
    Gompertz = safe_nls(
      Absorbancia ~ Bottom + A * exp(-exp(-k * (Tiempo - Ti))),
      start = list(Bottom = y_min, A = amp, k = 0.25, Ti = 14),
      lower = c(Bottom = 0, A = 0.1, k = 0.001, Ti = 0),
      upper = c(Bottom = 1.0, A = 3.0, k = 2.0, Ti = 56),
      model_name = "Gompertz"
    ),
    Weibull = safe_nls(
      Absorbancia ~ Bottom + A * (1 - exp(-(Tiempo / lambda)^beta)),
      start = list(Bottom = y_min, A = amp, lambda = 16, beta = 3),
      lower = c(Bottom = 0, A = 0.1, lambda = 1, beta = 0.2),
      upper = c(Bottom = 1.0, A = 3.0, lambda = 80, beta = 20),
      model_name = "Weibull"
    )
  )
}

nls_fits <- lapply(split(group_means, group_means$Grupo), fit_nls_models)
nls_summary <- bind_rows(lapply(names(nls_fits), function(grp) {
  bind_rows(lapply(nls_fits[[grp]], function(x) {
    if (is.null(x)) return(NULL)
    data.frame(Grupo = grp, Modelo = x$model, AIC = x$AIC, RSS = x$RSS, RMSE = x$RMSE, R2 = x$R2)
  }))
})) |>
  group_by(Grupo) |>
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE), best_model = delta_AIC == 0) |>
  ungroup() |>
  arrange(Grupo, AIC)
write_csv(nls_summary, file.path(tables_dir, "11_comparacion_modelos_no_lineales_promedio.csv"))

best_nls <- nls_summary |>
  filter(best_model) |>
  select(Grupo, Modelo)

extract_nls_params <- function(grp, model_name) {
  fit_obj <- nls_fits[[grp]][[model_name]]
  p <- fit_obj$params
  out <- data.frame(Grupo = grp, Modelo = model_name, parametro = names(p), valor = as.numeric(p))
  if (model_name == "Logistic") {
    derived <- data.frame(
      Grupo = grp, Modelo = model_name,
      parametro = c("amplitud", "max_velocity_abs_per_day"),
      valor = c(unname(p["Top"] - p["Bottom"]), unname((p["Top"] - p["Bottom"]) / (4 * p["s"])))
    )
    out <- bind_rows(out, derived)
  }
  if (model_name == "Gompertz") {
    derived <- data.frame(
      Grupo = grp, Modelo = model_name,
      parametro = c("Top_asymptote", "max_velocity_abs_per_day"),
      valor = c(unname(p["Bottom"] + p["A"]), unname(p["A"] * p["k"] / exp(1)))
    )
    out <- bind_rows(out, derived)
  }
  if (model_name == "Weibull") {
    derived <- data.frame(
      Grupo = grp, Modelo = model_name,
      parametro = c("Top_asymptote"),
      valor = c(unname(p["Bottom"] + p["A"]))
    )
    out <- bind_rows(out, derived)
  }
  out
}

nls_parameters <- bind_rows(mapply(
  extract_nls_params,
  as.character(best_nls$Grupo),
  as.character(best_nls$Modelo),
  SIMPLIFY = FALSE
))
write_csv(nls_parameters, file.path(tables_dir, "12_parametros_modelos_no_lineales.csv"))

x_pred <- seq(min(datos$Tiempo), max(datos$Tiempo), length.out = 300)
nls_predictions <- bind_rows(lapply(seq_len(nrow(best_nls)), function(i) {
  grp <- as.character(best_nls$Grupo[i])
  model_name <- as.character(best_nls$Modelo[i])
  fit <- nls_fits[[grp]][[model_name]]$fit
  nd <- data.frame(Tiempo = x_pred)
  data.frame(Grupo = grp, Modelo = model_name, Tiempo = x_pred, Absorbancia = as.numeric(predict(fit, nd)))
}))
nls_predictions$Grupo <- factor(nls_predictions$Grupo, levels = c("25ug", "50ug", "100ug", "Commercial"))
group_means$Grupo <- factor(group_means$Grupo, levels = c("25ug", "50ug", "100ug", "Commercial"))

pred_grid <- expand.grid(
  Tiempo = x_pred,
  Grupo = levels(datos$Grupo),
  KEEP.OUT.ATTRS = FALSE
) |>
  mutate(
    Animal = datos$Animal[1],
    TiempoF = factor(round(Tiempo / 7) * 7, levels = levels(datos$TiempoF)),
    TimeIndex = as.numeric(factor(round(Tiempo / 7) * 7, levels = tiempos))
  )

spline_plot_model <- fit_lme(Absorbancia ~ Grupo * ns(Tiempo, df = 4), datos, "REML")
pred_grid$fit <- predict(spline_plot_model, newdata = pred_grid, level = 0)

fig1 <- ggplot() +
  geom_line(
    data = datos,
    aes(Tiempo, Absorbancia, group = Animal, color = Grupo),
    linewidth = 0.25, alpha = 0.25
  ) +
  geom_point(
    data = datos,
    aes(Tiempo, Absorbancia, color = Grupo),
    size = 1.2, alpha = 0.42
  ) +
  geom_ribbon(
    data = mean_ci,
    aes(Tiempo, ymin = lower, ymax = upper, fill = Grupo),
    alpha = 0.13, linewidth = 0
  ) +
  geom_line(
    data = mean_ci,
    aes(Tiempo, mean, color = Grupo),
    linewidth = 0.75
  ) +
  geom_vline(xintercept = 21, linetype = "dashed", linewidth = 0.35, color = "grey40") +
  scale_color_manual(values = palette_nature, drop = FALSE) +
  scale_fill_manual(values = palette_nature, drop = FALSE) +
  scale_x_continuous(breaks = tiempos, expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(breaks = seq(0, 2.0, 0.5), expand = expansion(mult = c(0, 0.03))) +
  coord_cartesian(ylim = c(0, 2.15)) +
  labs(x = "Time after first immunization (days)", y = "ELISA absorbance (450 nm)") +
  theme_nature(8)
save_nature(fig1, "Figure_1_observed_ELISA_kinetics", 7.2, 3.6)

fig2_data <- animal_metrics |>
  filter(Grupo != "Control") |>
  select(Grupo, Animal, primary_peak_0_21, booster_delta_21_28, auc_total, persistence_ratio_56_peak) |>
  pivot_longer(
    cols = c(primary_peak_0_21, booster_delta_21_28, auc_total, persistence_ratio_56_peak),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  mutate(
    Metric = factor(
      Metric,
      levels = c("primary_peak_0_21", "booster_delta_21_28", "auc_total", "persistence_ratio_56_peak"),
      labels = c("Primary peak\n(days 0-21)", "Net rise\n(day 21-28)", "Total AUC\n(days 0-56)", "Persistence\n(day 56 / peak)")
    )
  )

fig2 <- ggplot(fig2_data, aes(Grupo, Value, color = Grupo, fill = Grupo)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.16, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.09, height = 0), size = 1.6, alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = palette_nature, drop = FALSE) +
  scale_fill_manual(values = palette_nature, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_nature(8) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
save_nature(fig2, "Figure_2_animal_level_kinetic_metrics", 7.2, 3.2)

fig3 <- ggplot() +
  geom_point(data = group_means, aes(Tiempo, Absorbancia, color = Grupo), size = 1.7) +
  geom_line(data = nls_predictions, aes(Tiempo, Absorbancia, color = Grupo), linewidth = 0.8) +
  facet_wrap(~ Grupo, nrow = 1) +
  scale_color_manual(values = palette_nature[-1], drop = FALSE) +
  scale_x_continuous(breaks = c(0, 14, 28, 42, 56), expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(breaks = seq(0, 2.0, 0.5), expand = expansion(mult = c(0, 0.03))) +
  coord_cartesian(ylim = c(0, 2.15)) +
  labs(x = "Time after first immunization (days)", y = "Mean absorbance (450 nm)") +
  theme_nature(8) +
  theme(legend.position = "none")
save_nature(fig3, "Figure_3_best_nonlinear_mean_fits", 7.2, 2.7)

sig_day <- pairwise_day |>
  filter(grepl("Commercial|25ug|50ug|100ug", contrast)) |>
  mutate(
    p_label = case_when(
      p.value < 0.001 ~ "<0.001",
      p.value < 0.01 ~ "<0.01",
      p.value < 0.05 ~ "<0.05",
      TRUE ~ "ns"
    ),
    neg_log10_p = -log10(pmax(p.value, 1e-4))
  )
write_csv(sig_day, file.path(tables_dir, "13_mapa_significancia_comparaciones_por_dia.csv"))

interpretation <- tibble::tribble(
  ~hallazgo, ~soporte_estadistico, ~implicacion_biologica,
  "Las vacunas recombinantes muestran seroconversion temprana fuerte antes del refuerzo.",
  "El modelo mixto por dia y las medias marginales separan los grupos recombinantes del control y de Commercial desde los dias 14-21.",
  "La formulacion recombinante induce una respuesta primaria rapida, compatible con reconocimiento temprano del antigeno vacunal.",
  "La vacuna Commercial presenta cinetica retardada en este diseno.",
  "Las curvas promedio no lineales estiman T50 mas tardio para Commercial que para 25ug, 50ug y 100ug.",
  "Su unica dosis produce una fase ascendente posterior al dia 21, por lo que no es biologicamente comparable a un esquema prime-boost en tiempos tempranos.",
  "La magnitud acumulada se evalua mejor con AUC por animal que con un unico punto final.",
  "ANOVA/Tukey de auc_total y auc_21_56 cuantifica diferencias integradas entre grupos vacunales.",
  "AUC resume cantidad y duracion de anticuerpos, una lectura mas cercana a exposicion humoral sostenida que el pico aislado.",
  "La persistencia al dia 56 es alta en los grupos recombinantes.",
  "La razon dia56/pico por animal captura mantenimiento de respuesta tras alcanzar meseta.",
  "Una razon cercana a 1 indica que la respuesta no decae marcadamente dentro de los 56 dias medidos.",
  "La inferencia para Commercial y Control debe interpretarse con cautela.",
  "Los tamanos muestrales son pequenos: Commercial n=3 y Control n=2.",
  "Los resultados son utiles para describir tendencias biologicas, pero conviene validarlos con mas animales si se busca una conclusion confirmatoria."
)
write_csv(interpretation, file.path(tables_dir, "14_implicacion_biologica_resumen.csv"))

sink(file.path(out_dir, "RESUMEN_ANALISIS_MEJORADO.txt"))
cat("Analisis mejorado de cinetica ELISA\n")
cat("====================================\n\n")
cat("Datos:\n")
print(table(datos$Grupo, datos$Tiempo))
cat("\nComparacion de modelos mixtos por AIC:\n")
print(model_comparison)
cat("\nANOVA del modelo mixto saturado por dia medido:\n")
print(anova(final_time_model))
cat("\nANOVA de metricas por animal:\n")
print(read_csv(file.path(tables_dir, "07_anova_metricas_por_animal.csv"), show_col_types = FALSE))
cat("\nTendencia lineal de dosis en recombinantes:\n")
print(dose_trends)
cat("\nMejores modelos no lineales por grupo:\n")
print(nls_summary |> filter(best_model))
cat("\nParametros no lineales derivados:\n")
print(nls_parameters)
cat("\nInterpretacion biologica:\n")
print(interpretation, n = Inf)
sink()

message("Analysis complete. Results written to: ", out_dir)



# 1. Asegurar que los factores incluyan al grupo Control en la base de la jerarquía
datos$Grupo <- factor(datos$Grupo, levels = c("Control", "25ug", "50ug", "100ug", "Commercial"))
mean_ci$Grupo <- factor(mean_ci$Grupo, levels = c("Control", "25ug", "50ug", "100ug", "Commercial"))
nls_predictions$Grupo <- factor(nls_predictions$Grupo, levels = c("25ug", "50ug", "100ug", "Commercial"))

fig3_final_con_control <- ggplot() +
  # A. PUNTOS INDIVIDUALES (Fondo): Todos los animales, incluido el Control, con sutil jitter
  geom_jitter(
    data = datos,
    aes(x = Tiempo, y = Absorbancia, color = Grupo),
    size = 0.9,
    alpha = 0.15,
    width = 0.6,
    height = 0
  ) +
  
  # B. LÍNEA EMPÍRICA DEL GRUPO CONTROL: Conecta los promedios observados del Control punto a punto
  geom_line(
    data = mean_ci |> filter(Grupo == "Control"),
    aes(x = Tiempo, y = mean, color = Grupo),
    linewidth = 0.6,
    linetype = "twodash", # Estilo de línea diferente para dejar en claro que no es un modelo NLS
    alpha = 0.7
  ) +
  
  # C. LÍNEAS MODELADAS NO LINEALES (Frente): Curvas continuas para los grupos vacunados
  geom_line(
    data = nls_predictions, 
    aes(x = Tiempo, y = Absorbancia, color = Grupo), 
    linewidth = 0.95, 
    alpha = 0.95
  ) +
  
  # D. PUNTOS PROMEDIO (Medio): Marcadores de tendencia central para TODOS los grupos
  geom_point(
    data = mean_ci, 
    aes(x = Tiempo, y = mean, color = Grupo), 
    size = 2.0, 
    alpha = 0.75, 
    shape = 16
  ) +
  
  # E. LÍNEA VERTICAL DEL BOOSTER (Día 21)
  geom_vline(
    xintercept = 21, 
    linetype = "dashed", 
    linewidth = 0.35, 
    color = "grey40",
    alpha = 0.5
  ) +
  
  # F. ANOTACIÓN TEXTO "Booster"
  annotate(
    "text", 
    x = 22, 
    y = 2.0, 
    label = "Booster", 
    hjust = 0, 
    vjust = 1, 
    size = 2.5, 
    fontface = "italic", 
    color = "grey30"
  ) +
  
  # G. PALETAS Y ESCALAS: Usamos 'palette_nature' completa (con Control incluido)
  scale_color_manual(values = palette_nature) + 
  scale_x_continuous(
    breaks = c(0, 7, 14, 21, 28, 35, 42, 49, 56), 
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    breaks = seq(0, 2.0, 0.5), 
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_cartesian(ylim = c(0, 2.15)) +
  
  labs(
    x = "Time after first immunization (days)", 
    y = "Anti-E2 IgG (Absorbance, 450 nm)",
    color = "Vaccine Group"
  ) +
  
  theme_nature(8) +
  theme(
    legend.position = "top",
    legend.justification = "center"
  )

# Guardar la versión definitiva apta para publicación
save_nature(fig3_final_con_control, "Figure_3_unified_nonlinear_fits_with_control", 5.4, 4.0)






library(ggplot2)
library(dplyr)
library(tidyr)
library(emmeans)
library(purrr)
library(patchwork)

# 1. Asegurar el orden de los factores
animal_metrics$Grupo <- factor(animal_metrics$Grupo, levels = c("Control", "25ug", "50ug", "100ug", "Commercial"))

# Definir las 4 métricas con sus respectivos labels de eje Y específicos
metrics_to_plot <- c("primary_peak_0_21", "booster_delta_21_28", "auc_total", "persistence_ratio_56_peak")

# Modificación en la lista de metadatos de los ejes Y
metric_metadata <- list(
  "primary_peak_0_21" = list(
    title = "Primary peak\n(days 0-21)",
    ylabel = "Anti-E2 Antibodies (Absorbance, 450 nm)"
  ),
  "booster_delta_21_28" = list(
    title = "Net rise\n(day 21-28)",
    ylabel = "Net Increase (Abs D28 - Abs D21)" # Fórmula añadida
  ),
  "auc_total" = list(
    title = "Total AUC\n(days 0-56)",
    ylabel = "Total Antibody Exposure (Absorbance x Days)"
  ),
  "persistence_ratio_56_peak" = list(
    title = "Persistence\n(day 56 / peak)",
    ylabel = "Persistence Ratio (Abs D56 / Max Peak)" # Fórmula añadida
  )
)

# [El resto del código de la función 'generate_individual_boxplot' y el ensamblado final se mantiene exactamente igual]

# 2. Función para generar cada gráfico con su propia unidad en el eje Y
generate_individual_boxplot <- function(metric_name, meta) {
  
  # Filtrar datos crudos para este panel
  plot_data <- animal_metrics |> 
    select(Grupo, Value = .data[[metric_name]]) |> 
    filter(!is.na(Value))
  
  # Calcular estadística EXCLUYENDO al grupo Control
  stats_data <- plot_data |> filter(Grupo != "Control") |> droplevels()
  fit <- lm(Value ~ Grupo, data = stats_data)
  emm <- emmeans(fit, ~ Grupo)
  pairs_df <- as.data.frame(summary(pairs(emm, adjust = "tukey")))
  
  max_y <- max(plot_data$Value, na.rm = TRUE)
  range_y <- max_y - min(plot_data$Value, na.rm = TRUE)
  
  # Configurar barras de comparación vs Commercial (X: Control=1, 25ug=2, 50ug=3, 100ug=4, Commercial=5)
  comp_metadata <- list(
    list(g1 = "25ug",  x1 = 2, x2 = 5, row = 1, pair_name = "25ug - Commercial"),
    list(g1 = "50ug",  x1 = 3, x2 = 5, row = 2, pair_name = "50ug - Commercial"),
    list(g1 = "100ug", x1 = 4, x2 = 5, row = 3, pair_name = "100ug - Commercial")
  )
  
  lines_and_labels <- map_df(comp_metadata, function(m) {
    p_val <- pairs_df$p.value[pairs_df$contrast == m$pair_name | pairs_df$contrast == paste0("Commercial - ", m$g1)]
    if(length(p_val) == 0) {
      p_val <- pairs_df |> filter(grepl(m$g1, contrast) & grepl("Commercial", contrast)) |> pull(p.value)
    }
    
    p_formatted <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", p_val))
    y_bar <- max_y + (range_y * 0.09 * m$row)
    
    data.frame(
      x = m$x1, xend = m$x2, y = y_bar,
      x_text = (m$x1 + m$x2) / 2, y_text = y_bar + (range_y * 0.03),
      p_label = p_formatted
    )
  })
  
  ylim_max <- max(lines_and_labels$y_text) + (range_y * 0.05)
  
  # 3. Construir el gráfico con su label de eje Y personalizado
  p <- ggplot(plot_data, aes(x = Grupo, y = Value, color = Grupo)) +
    geom_boxplot(aes(fill = Grupo), width = 0.5, outlier.shape = NA, alpha = 0.14, linewidth = 0.4) +
    geom_point(position = position_jitter(width = 0.1, height = 0), size = 1.3, alpha = 0.75) +
    
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
    coord_cartesian(ylim = c(min(plot_data$Value) - (range_y * 0.05), ylim_max)) +
    
    labs(
      title = meta$title,
      x = NULL, 
      y = meta$ylabel # Aquí se aplica el eje Y correcto para cada panel
    ) +
    
    theme_nature(8) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.title.y = element_text(size = 7.5, face = "plain"),
      plot.margin = margin(t = 5, r = 8, b = 5, l = 5)
    )
  
  return(p)
}

# 4. Generar los paneles independientes mapeando la lista de metadatos
plot_list <- imap(metric_metadata, ~ generate_individual_boxplot(.y, .x))

# 5. Combinar manteniendo la total independencia de sus ejes Y
fig2_final_unidades_ok <- patchwork::wrap_plots(plot_list, nrow = 1) + 
  patchwork::plot_annotation(
    theme = theme(plot.margin = margin(10, 10, 10, 10))
  )

# 6. Guardar la figura corregida
save_nature(fig2_final_unidades_ok, "Figure_2_animal_metrics_correct_units", 9.2, 3.8)

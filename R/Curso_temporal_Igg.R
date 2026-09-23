# ============================================================
# ANALISIS TEMPORAL DE IgG ANTI-E2
# Proyecto: Manuscrito-BVDV-challenge
# Script: 01_IgG_temporal.R
#
# Entrada:
#   data/processed/Curso_temporal_IgG.xlsx
#
# Salidas:
#   results/tables/Curso_temporal_IgG_resultados.xlsx
#   results/figures/
#   results/models/
#
# ============================================================


# ============================================================
# 1. PAQUETES
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "openxlsx",
  "nlme",
  "emmeans",
  "broom",
  "patchwork",
  "cowplot",
  "pracma",
  "scales",
  "splines",
  "purrr"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(openxlsx)
  library(nlme)
  library(emmeans)
  library(broom)
  library(patchwork)
  library(cowplot)
  library(pracma)
  library(scales)
  library(splines)
  library(purrr)
})


# ============================================================
# 2. SEMILLA
# ============================================================

set.seed(20260620)


# ============================================================
# 3. IDENTIFICAR DIRECTORIO DEL PROYECTO
# ============================================================

find_project_root <- function() {
  
  current <- normalizePath(getwd(), mustWork = TRUE)
  
  candidates <- unique(c(
    current,
    dirname(current),
    dirname(dirname(current))
  ))
  
  for (path in candidates) {
    
    if (file.exists(
      file.path(path, "Manuscrito-BVDV-challenge.Rproj")
    )) {
      return(path)
    }
  }
  
  stop(
    paste0(
      "No se pudo encontrar el directorio raíz del proyecto.\n",
      "Asegúrate de ejecutar este script desde el proyecto ",
      "'Manuscrito-BVDV-challenge'."
    )
  )
}

project_dir <- find_project_root()

cat("Directorio del proyecto:\n")
cat(project_dir, "\n\n")


# ============================================================
# 4. DIRECTORIOS
# ============================================================

data_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

results_dir <- file.path(
  project_dir,
  "results"
)

tables_dir <- file.path(
  results_dir,
  "tables"
)

figures_dir <- file.path(
  results_dir,
  "figures"
)

models_dir <- file.path(
  results_dir,
  "models"
)

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 5. ARCHIVO DE ENTRADA
# ============================================================

input_file <- file.path(
  data_dir,
  "Curso_temporal_IgG.xlsx"
)

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "\nNo se encontró el archivo:\n",
      input_file,
      "\n\n",
      "Coloca 'Curso_temporal_IgG.xlsx' en:\n",
      data_dir
    )
  )
}


# ============================================================
# 6. LEER DATOS
# ============================================================

cat("Hojas disponibles en el Excel:\n")

print(excel_sheets(input_file))

raw_data <- read_excel(
  input_file,
  sheet = 1
)

cat("\nColumnas encontradas:\n")
print(names(raw_data))


# ============================================================
# 7. VALIDAR COLUMNAS
# ============================================================

required_columns <- c(
  "Animal",
  "Grupo",
  "Tiempo",
  "Absorbancia"
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "\nFaltan las siguientes columnas obligatorias:\n",
      paste(missing_columns, collapse = ", "),
      "\n\nLas columnas esperadas son:\n",
      paste(required_columns, collapse = ", ")
    )
  )
}


# ============================================================
# 8. PREPARAR DATOS
# ============================================================

datos <- raw_data %>%
  
  mutate(
    
    Animal = factor(Animal),
    
    Grupo = as.character(Grupo),
    
    # Estandarizar nombres
    Grupo = case_when(
      Grupo == "Comercial" ~ "Commercial",
      Grupo == "Commercial" ~ "Commercial",
      TRUE ~ Grupo
    ),
    
    Grupo = factor(
      Grupo,
      levels = c(
        "Control",
        "25ug",
        "50ug",
        "100ug",
        "Commercial"
      )
    ),
    
    Tiempo = as.numeric(Tiempo),
    
    Absorbancia = as.numeric(Absorbancia)
  ) %>%
  
  filter(
    !is.na(Animal),
    !is.na(Grupo),
    !is.na(Tiempo),
    !is.na(Absorbancia)
  ) %>%
  
  mutate(
    
    TiempoF = factor(
      Tiempo,
      levels = sort(unique(Tiempo))
    ),
    
    TimeIndex = as.numeric(
      factor(
        Tiempo,
        levels = sort(unique(Tiempo))
      )
    )
  )


# ============================================================
# 9. COMPROBACIÓN DE DATOS
# ============================================================

cat("\nNúmero de observaciones:", nrow(datos), "\n")

cat("\nNúmero de animales:",
    n_distinct(datos$Animal),
    "\n")

cat("\nGrupos:\n")
print(table(datos$Grupo))

cat("\nTiempos:\n")
print(sort(unique(datos$Tiempo)))


# ============================================================
# 10. DATOS VACUNADOS
# ============================================================

datos_vac <- datos %>%
  filter(
    Grupo != "Control"
  ) %>%
  droplevels()

tiempos <- sort(
  unique(datos$Tiempo)
)


# ============================================================
# 11. PALETA
# ============================================================

palette_nature <- c(
  "Control" = "#4D4D4D",
  "25ug" = "#004B87",
  "50ug" = "#E69F00",
  "100ug" = "#2ECC71",
  "Commercial" = "#D55E00"
)


# ============================================================
# 12. TEMA
# ============================================================

theme_nature <- function(base_size = 9) {
  
  theme_classic(
    base_size = base_size
  ) +
    
    theme(
      
      plot.title = element_text(
        face = "bold",
        size = base_size + 1
      ),
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text = element_text(
        color = "black"
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      legend.position = "right",
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
}


# ============================================================
# 13. FUNCIÓN PARA GUARDAR FIGURAS
# ============================================================

save_nature <- function(
    plot,
    filename,
    width = 6,
    height = 4
) {
  
  ggsave(
    file.path(
      figures_dir,
      paste0(filename, ".pdf")
    ),
    plot,
    width = width,
    height = height,
    units = "in"
  )
  
  ggsave(
    file.path(
      figures_dir,
      paste0(filename, ".tiff")
    ),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    file.path(
      figures_dir,
      paste0(filename, ".png")
    ),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
}


# ============================================================
# 14. RESUMEN POR GRUPO Y TIEMPO
# ============================================================

mean_ci <- datos %>%
  
  group_by(
    Grupo,
    Tiempo
  ) %>%
  
  summarise(
    
    n = n(),
    
    mean = mean(
      Absorbancia,
      na.rm = TRUE
    ),
    
    sd = sd(
      Absorbancia,
      na.rm = TRUE
    ),
    
    se = sd / sqrt(n),
    
    ci95 = qt(
      0.975,
      df = pmax(n - 1, 1)
    ) * se,
    
    lower = mean - ci95,
    
    upper = mean + ci95,
    
    .groups = "drop"
  )


# ============================================================
# 15. MÉTRICAS INDIVIDUALES
# ============================================================

animal_metrics <- datos_vac %>%
  
  group_by(
    Animal,
    Grupo
  ) %>%
  
  summarise(
    
    baseline = ifelse(
      any(Tiempo == 0),
      Absorbancia[Tiempo == 0][1],
      NA_real_
    ),
    
    day14 = ifelse(
      any(Tiempo == 14),
      Absorbancia[Tiempo == 14][1],
      NA_real_
    ),
    
    day21 = ifelse(
      any(Tiempo == 21),
      Absorbancia[Tiempo == 21][1],
      NA_real_
    ),
    
    day28 = ifelse(
      any(Tiempo == 28),
      Absorbancia[Tiempo == 28][1],
      NA_real_
    ),
    
    day56 = ifelse(
      any(Tiempo == 56),
      Absorbancia[Tiempo == 56][1],
      NA_real_
    ),
    
    peak = max(
      Absorbancia,
      na.rm = TRUE
    ),
    
    peak_day = Tiempo[
      which.max(Absorbancia)
    ],
    
    delta_D28_D21 = day28 - day21,
    
    persistence = ifelse(
      !is.na(day56) & peak > 0,
      day56 / peak,
      NA_real_
    ),
    
    AUC = if (
      sum(!is.na(Absorbancia)) >= 2
    ) {
      
      pracma::trapz(
        Tiempo,
        Absorbancia
      )
      
    } else {
      
      NA_real_
    },
    
    .groups = "drop"
  )


# ============================================================
# 16. RESUMEN DE MÉTRICAS POR GRUPO
# ============================================================

metric_summary <- animal_metrics %>%
  
  pivot_longer(
    
    cols = c(
      peak,
      delta_D28_D21,
      AUC,
      persistence
    ),
    
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  
  group_by(
    Grupo,
    Metric
  ) %>%
  
  summarise(
    
    n = sum(!is.na(Value)),
    
    mean = mean(
      Value,
      na.rm = TRUE
    ),
    
    sd = sd(
      Value,
      na.rm = TRUE
    ),
    
    se = sd / sqrt(n),
    
    lower = mean -
      qt(
        0.975,
        df = pmax(n - 1, 1)
      ) * se,
    
    upper = mean +
      qt(
        0.975,
        df = pmax(n - 1, 1)
      ) * se,
    
    .groups = "drop"
  )


# ============================================================
# 17. MODELOS MIXTOS
# ============================================================

fit_lme <- function(
    expr,
    data,
    method = "REML",
    correlation = NULL
) {
  
  tryCatch({
    
    model <- nlme::lme(
      fixed = expr,
      random = ~1 | Animal,
      data = data,
      method = method,
      correlation = correlation,
      na.action = na.omit,
      keep.data = TRUE
    )
    
    return(model)
    
  }, error = function(e) {
    
    message(
      "No se pudo ajustar el modelo: ",
      conditionMessage(e)
    )
    
    return(NULL)
  })
}


# ============================================================
# 18. MODELOS CANDIDATOS
# ============================================================

model_candidates <- list(
  
  measured_time = fit_lme(
    Absorbancia ~ Grupo * TiempoF,
    datos,
    method = "ML"
  ),
  
  measured_time_AR1 = fit_lme(
    Absorbancia ~ Grupo * TiempoF,
    datos,
    method = "ML",
    correlation = corAR1(
      form = ~TimeIndex | Animal
    )
  ),
  
  spline_df3 = fit_lme(
    Absorbancia ~ Grupo * ns(
      Tiempo,
      df = 3
    ),
    datos,
    method = "ML"
  ),
  
  spline_df4 = fit_lme(
    Absorbancia ~ Grupo * ns(
      Tiempo,
      df = 4
    ),
    datos,
    method = "ML"
  ),
  
  cubic = fit_lme(
    Absorbancia ~ Grupo * poly(
      Tiempo,
      3
    ),
    datos,
    method = "ML"
  )
)


# ============================================================
# 19. COMPARACIÓN DE MODELOS
# ============================================================

model_comparison <- map_dfr(
  names(model_candidates),
  function(name) {
    
    model <- model_candidates[[name]]
    
    if (is.null(model)) {
      
      return(
        tibble(
          Model = name,
          AIC = NA_real_,
          BIC = NA_real_,
          logLik = NA_real_
        )
      )
    }
    
    tibble(
      
      Model = name,
      
      AIC = AIC(model),
      
      BIC = BIC(model),
      
      logLik = as.numeric(
        logLik(model)
      )
    )
  }
) %>%
  
  arrange(AIC)


# ============================================================
# 20. MODELO FINAL
# ============================================================

final_time_model <- fit_lme(
  Absorbancia ~ Grupo * TiempoF,
  datos,
  method = "REML"
)

if (!is.null(final_time_model)) {
  
  saveRDS(
    final_time_model,
    file.path(
      models_dir,
      "IgG_modelo_mixto_tiempo_REML.rds"
    )
  )
}


# ============================================================
# 21. ANOVA DEL MODELO FINAL
# ============================================================

anova_final <- if (
  !is.null(final_time_model)
) {
  
  as.data.frame(
    anova(final_time_model)
  ) %>%
    tibble::rownames_to_column(
      "Effect"
    )
  
} else {
  
  tibble()
}


# ============================================================
# 22. EMMEANS POR TIEMPO
# ============================================================

emmeans_time <- if (
  !is.null(final_time_model)
) {
  
  emmeans(
    final_time_model,
    ~ Grupo | TiempoF,
    data = datos
  ) %>%
    as.data.frame()
  
} else {
  
  tibble()
}


# ============================================================
# 23. COMPARACIONES ENTRE GRUPOS
# ============================================================

pairwise_time <- if (
  !is.null(final_time_model)
) {
  
  emmeans(
    final_time_model,
    ~ Grupo | TiempoF,
    data = datos
  ) %>%
    pairs(
      adjust = "tukey"
    ) %>%
    as.data.frame()
  
} else {
  
  tibble()
}

# ============================================================
# 24. MÉTRICAS: MODELOS POR VARIABLE
# ============================================================

metric_names <- c(
  "peak",
  "delta_D28_D21",
  "AUC",
  "persistence"
)


metric_models <- list()


metric_tests <- map_dfr(
  metric_names,
  function(metric) {
    
    formula <- as.formula(
      paste(
        metric,
        "~ Grupo"
      )
    )
    
    data_metric <- animal_metrics %>%
      
      filter(
        !is.na(
          .data[[metric]]
        )
      )
    
    if (
      nrow(data_metric) < 3
    ) {
      return(
        tibble()
      )
    }
    
    model <- lm(
      formula,
      data = data_metric
    )
    
    metric_models[[metric]] <<- model
    
    broom::tidy(
      anova(model)
    ) %>%
      
      mutate(
        Metric = metric,
        .before = 1
      )
  }
)


# ============================================================
# 25. EMMEANS PARA MÉTRICAS
# ============================================================

metric_emmeans <- map_dfr(
  metric_names,
  function(metric) {
    
    data_metric <- animal_metrics %>%
      
      filter(
        !is.na(
          .data[[metric]]
        )
      )
    
    if (
      nrow(data_metric) < 3
    ) {
      return(
        tibble()
      )
    }
    
    formula <- as.formula(
      paste(
        metric,
        "~ Grupo"
      )
    )
    
    model <- lm(
      formula,
      data = data_metric
    )
    
    emmeans(
      model,
      ~Grupo
    ) %>%
      
      pairs(
        adjust = "tukey"
      ) %>%
      
      as.data.frame() %>%
      
      mutate(
        Metric = metric,
        .before = 1
      )
  }
)


# ============================================================
# 26. TENDENCIA DE DOSIS
# ============================================================

dose_data <- animal_metrics %>%
  
  filter(
    Grupo %in% c(
      "25ug",
      "50ug",
      "100ug"
    )
  ) %>%
  
  mutate(
    
    Dose = case_when(
      
      Grupo == "25ug" ~ 25,
      
      Grupo == "50ug" ~ 50,
      
      Grupo == "100ug" ~ 100
      
    )
  )


dose_tests <- map_dfr(
  metric_names,
  function(metric) {
    
    d <- dose_data %>%
      
      filter(
        !is.na(
          .data[[metric]]
        )
      )
    
    if (
      nrow(d) < 3
    ) {
      return(
        tibble()
      )
    }
    
    formula <- as.formula(
      paste(
        metric,
        "~ Dose"
      )
    )
    
    model <- lm(
      formula,
      data = d
    )
    
    broom::tidy(
      model
    ) %>%
      
      mutate(
        Metric = metric,
        .before = 1
      )
  }
)


# ============================================================
# 27. MODELOS NO LINEALES SOBRE LAS MEDIAS
# ============================================================

nls_data <- mean_ci %>%
  
  filter(
    !is.na(mean)
  )


# ------------------------------------------------------------
# Modelo logístico
# ------------------------------------------------------------

logistic_model <- tryCatch(
  
  nls(
    mean ~
      A /
      (
        1 +
          exp(
            -(Tiempo - T50) / k
          )
      ),
    
    data = nls_data,
    
    start = list(
      A = max(
        nls_data$mean,
        na.rm = TRUE
      ),
      T50 = median(
        nls_data$Tiempo
      ),
      k = 10
    ),
    
    control = nls.control(
      maxiter = 1000
    )
  ),
  
  error = function(e) NULL
)


# ------------------------------------------------------------
# Modelo Gompertz
# ------------------------------------------------------------

gompertz_model <- tryCatch(
  
  nls(
    mean ~
      A *
      exp(
        -exp(
          -(Tiempo - T50) / k
        )
      ),
    
    data = nls_data,
    
    start = list(
      A = max(
        nls_data$mean,
        na.rm = TRUE
      ),
      T50 = median(
        nls_data$Tiempo
      ),
      k = 10
    ),
    
    control = nls.control(
      maxiter = 1000
    )
  ),
  
  error = function(e) NULL
)


# ------------------------------------------------------------
# Modelo Weibull
# ------------------------------------------------------------

weibull_model <- tryCatch(
  
  nls(
    mean ~
      A *
      (
        1 -
          exp(
            -(Tiempo / lambda)^k
          )
      ),
    
    data = nls_data,
    
    start = list(
      A = max(
        nls_data$mean,
        na.rm = TRUE
      ),
      lambda = median(
        nls_data$Tiempo[
          nls_data$Tiempo > 0
        ]
      ),
      k = 2
    ),
    
    control = nls.control(
      maxiter = 1000
    )
  ),
  
  error = function(e) NULL
)


# ============================================================
# 28. COMPARACIÓN DE MODELOS NLS
# ============================================================

nls_models <- list(
  
  Logistic = logistic_model,
  
  Gompertz = gompertz_model,
  
  Weibull = weibull_model
)


nls_comparison <- map_dfr(
  names(nls_models),
  function(name) {
    
    model <- nls_models[[name]]
    
    if (is.null(model)) {
      
      return(
        tibble(
          Model = name,
          AIC = NA_real_,
          RSS = NA_real_,
          RMSE = NA_real_,
          R2_pseudo = NA_real_
        )
      )
    }
    
    residuals_model <- residuals(model)
    
    fitted_model <- fitted(model)
    
    rss <- sum(
      residuals_model^2
    )
    
    rmse <- sqrt(
      mean(
        residuals_model^2
      )
    )
    
    tss <- sum(
      (
        nls_data$mean -
          mean(nls_data$mean)
      )^2
    )
    
    r2 <- 1 - rss / tss
    
    tibble(
      
      Model = name,
      
      AIC = AIC(model),
      
      RSS = rss,
      
      RMSE = rmse,
      
      R2_pseudo = r2
    )
  }
) %>%
  
  arrange(AIC)


# ============================================================
# 29. PARÁMETROS NLS
# ============================================================

nls_parameters <- map_dfr(
  names(nls_models),
  function(name) {
    
    model <- nls_models[[name]]
    
    if (is.null(model)) {
      
      return(
        tibble()
      )
    }
    
    broom::tidy(model) %>%
      
      mutate(
        Model = name,
        .before = 1
      )
  }
)


# ============================================================
# 30. PREDICCIONES NLS
# ============================================================

prediction_grid <- tibble(
  
  Tiempo = seq(
    min(nls_data$Tiempo),
    max(nls_data$Tiempo),
    length.out = 200
  )
)


nls_predictions <- map_dfr(
  names(nls_models),
  function(name) {
    
    model <- nls_models[[name]]
    
    if (is.null(model)) {
      return(tibble())
    }
    
    prediction_grid %>%
      
      mutate(
        
        Predicted = predict(
          model,
          newdata = prediction_grid
        ),
        
        Model = name
      )
  }
)


# ============================================================
# 31. FIGURA 1
# CINÉTICA OBSERVADA
# ============================================================

fig1 <- ggplot(
  
  mean_ci,
  
  aes(
    x = Tiempo,
    y = mean,
    color = Grupo,
    group = Grupo
  )
  
) +
  
  geom_ribbon(
    
    aes(
      ymin = lower,
      ymax = upper,
      fill = Grupo
    ),
    
    alpha = 0.15,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 1.8
  ) +
  
  scale_color_manual(
    values = palette_nature,
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = palette_nature,
    drop = FALSE
  ) +
  
  labs(
    
    title = "Temporal kinetics of anti-E2 IgG",
    
    x = "Days post-immunization",
    
    y = "Absorbance"
  ) +
  
  theme_nature()


save_nature(
  fig1,
  "IgG_Figure_1_temporal_kinetics",
  6,
  4
)


# ============================================================
# 32. FIGURA 2
# MÉTRICAS INDIVIDUALES
# ============================================================

metric_labels <- c(
  
  peak =
    "Primary peak\nAnti-E2 IgG absorbance",
  
  delta_D28_D21 =
    "Net rise\nD28 - D21",
  
  AUC =
    "Total AUC\nAbsorbance × days",
  
  persistence =
    "Persistence\nD56 / peak"
)


generate_metric_plot <- function(
    metric_name
) {
  
  d <- animal_metrics %>%
    
    filter(
      !is.na(
        .data[[metric_name]]
      )
    )
  
  ggplot(
    d,
    aes(
      x = Grupo,
      y = .data[[metric_name]],
      fill = Grupo
    )
  ) +
    
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      width = 0.08,
      size = 1.6
    ) +
    
    scale_fill_manual(
      values = palette_nature,
      drop = FALSE
    ) +
    
    labs(
      
      x = NULL,
      
      y = metric_labels[
        metric_name
      ]
    ) +
    
    theme_nature() +
    
    theme(
      legend.position = "none"
    )
}


fig2_peak <- generate_metric_plot(
  "peak"
)

fig2_delta <- generate_metric_plot(
  "delta_D28_D21"
)

fig2_auc <- generate_metric_plot(
  "AUC"
)

fig2_persistence <- generate_metric_plot(
  "persistence"
)


fig2 <- (
  
  fig2_peak |
    
    fig2_delta |
    
    fig2_auc |
    
    fig2_persistence
  
) +
  
  plot_annotation(
    title =
      "Individual-level kinetic metrics"
  )


save_nature(
  fig2,
  "IgG_Figure_2_individual_metrics",
  10,
  3.2
)


# ============================================================
# 33. FIGURA 3
# MODELOS NO LINEALES
# ============================================================

fig3_data <- mean_ci %>%
  
  filter(
    !is.na(mean)
  )


fig3 <- ggplot() +
  
  geom_point(
    
    data = fig3_data,
    
    aes(
      x = Tiempo,
      y = mean,
      color = Grupo
    ),
    
    size = 1.5
  ) +
  
  geom_line(
    
    data = fig3_data,
    
    aes(
      x = Tiempo,
      y = mean,
      color = Grupo,
      group = Grupo
    ),
    
    alpha = 0.4
  ) +
  
  geom_line(
    
    data = nls_predictions,
    
    aes(
      x = Tiempo,
      y = Predicted,
      linetype = Model
    ),
    
    linewidth = 0.8
  ) +
  
  scale_color_manual(
    values = palette_nature,
    drop = FALSE
  ) +
  
  labs(
    
    title =
      "Nonlinear models of temporal IgG kinetics",
    
    x =
      "Days post-immunization",
    
    y =
      "Mean absorbance",
    
    linetype =
      "Model"
  ) +
  
  theme_nature()


save_nature(
  fig3,
  "IgG_Figure_3_nonlinear_models",
  6,
  4
)


# ============================================================
# 34. TABLA DE DATOS
# ============================================================

datos_export <- datos %>%
  
  mutate(
    Grupo = as.character(Grupo)
  )


# ============================================================
# 35. TABLA README DEL EXCEL
# ============================================================

README_excel <- tibble(
  
  Sheet = c(
    "README",
    "Datos",
    "Resumen",
    "Metricas_animal",
    "Metricas_resumen",
    "Modelo_mixto",
    "EMMeans",
    "Comparaciones_tiempo",
    "Metricas_ANOVA",
    "Metricas_comparaciones",
    "Dosis",
    "NLS_comparacion",
    "NLS_parametros",
    "NLS_predicciones"
  ),
  
  Description = c(
    
    "Descripción de las hojas del archivo",
    
    "Datos procesados utilizados para el análisis",
    
    "Media, SD, SE e IC95% por grupo y tiempo",
    
    "Métricas calculadas para cada animal",
    
    "Resumen de métricas por grupo",
    
    "Comparación de modelos mixtos longitudinales",
    
    "Estimated marginal means del modelo final",
    
    "Comparaciones Tukey entre grupos por tiempo",
    
    "ANOVA de las métricas individuales",
    
    "Comparaciones entre grupos para las métricas",
    
    "Análisis de tendencia según dosis de antígeno",
    
    "Comparación de modelos Logistic, Gompertz y Weibull",
    
    "Parámetros estimados de los modelos NLS",
    
    "Predicciones de los modelos no lineales"
  )
)


# ============================================================
# 36. CREAR ARCHIVO EXCEL DE RESULTADOS
# ============================================================

output_excel <- file.path(
  
  tables_dir,
  
  "Curso_temporal_IgG_resultados.xlsx"
)


# Eliminar archivo anterior si existe

if (file.exists(output_excel)) {
  
  file.remove(output_excel)
}


wb <- createWorkbook()


# ------------------------------------------------------------
# Función para agregar hojas
# ------------------------------------------------------------

add_sheet <- function(
    wb,
    sheet_name,
    data
) {
  
  addWorksheet(
    wb,
    sheetName = sheet_name
  )
  
  if (nrow(data) > 0) {
    
    writeData(
      wb,
      sheet = sheet_name,
      x = data,
      withFilter = TRUE
    )
    
  } else {
    
    writeData(
      wb,
      sheet = sheet_name,
      x = data
    )
  }
  
  freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:ncol(data),
    widths = "auto"
  )
}


# ============================================================
# 37. AGREGAR TODAS LAS HOJAS
# ============================================================

add_sheet(
  wb,
  "README",
  README_excel
)

add_sheet(
  wb,
  "Datos",
  datos_export
)

add_sheet(
  wb,
  "Resumen",
  mean_ci
)

add_sheet(
  wb,
  "Metricas_animal",
  animal_metrics
)

add_sheet(
  wb,
  "Metricas_resumen",
  metric_summary
)

add_sheet(
  wb,
  "Modelo_mixto",
  model_comparison
)

add_sheet(
  wb,
  "EMMeans",
  emmeans_time
)

add_sheet(
  wb,
  "Comparaciones_tiempo",
  pairwise_time
)

add_sheet(
  wb,
  "Metricas_ANOVA",
  metric_tests
)

add_sheet(
  wb,
  "Metricas_comparaciones",
  metric_emmeans
)

add_sheet(
  wb,
  "Dosis",
  dose_tests
)

add_sheet(
  wb,
  "NLS_comparacion",
  nls_comparison
)

add_sheet(
  wb,
  "NLS_parametros",
  nls_parameters
)

add_sheet(
  wb,
  "NLS_predicciones",
  nls_predictions
)


# ============================================================
# 38. ESTILO DEL EXCEL
# ============================================================

header_style <- createStyle(
  
  textDecoration = "bold",
  
  halign = "center",
  
  border = "Bottom"
)


for (
  sheet in names(wb)
) {
  
  addStyle(
    
    wb,
    
    sheet = sheet,
    
    style = header_style,
    
    rows = 1,
    
    cols = 1:max(
      1,
      ncol(
        readWorkbook(
          wb,
          sheet = sheet,
          rows = 1
        )
      )
    ),
    
    gridExpand = TRUE
  )
}


# ============================================================
# 39. GUARDAR EXCEL
# ============================================================

saveWorkbook(
  
  wb,
  
  output_excel,
  
  overwrite = TRUE
)


# ============================================================
# 40. GUARDAR MODELOS NLS
# ============================================================

if (!is.null(logistic_model)) {
  
  saveRDS(
    
    logistic_model,
    
    file.path(
      models_dir,
      "IgG_NLS_Logistic.rds"
    )
  )
}


if (!is.null(gompertz_model)) {
  
  saveRDS(
    
    gompertz_model,
    
    file.path(
      models_dir,
      "IgG_NLS_Gompertz.rds"
    )
  )
}


if (!is.null(weibull_model)) {
  
  saveRDS(
    
    weibull_model,
    
    file.path(
      models_dir,
      "IgG_NLS_Weibull.rds"
    )
  )
}


# ============================================================
# 41. GUARDAR RESUMEN DE EJECUCIÓN
# ============================================================

analysis_info <- tibble(
  
  Item = c(
    
    "Analysis",
    
    "Input file",
    
    "Number of observations",
    
    "Number of animals",
    
    "Groups",
    
    "Time points",
    
    "Output Excel",
    
    "Figures directory",
    
    "Models directory",
    
    "Analysis date"
  ),
  
  Value = c(
    
    "Temporal anti-E2 IgG analysis",
    
    input_file,
    
    nrow(datos),
    
    n_distinct(datos$Animal),
    
    paste(
      levels(datos$Grupo),
      collapse = ", "
    ),
    
    paste(
      sort(unique(datos$Tiempo)),
      collapse = ", "
    ),
    
    output_excel,
    
    figures_dir,
    
    models_dir,
    
    as.character(
      Sys.Date()
    )
  )
)


write_csv(
  
  analysis_info,
  
  file.path(
    tables_dir,
    "IgG_analysis_info.csv"
  )
)


# ============================================================
# 42. MENSAJE FINAL
# ============================================================

cat("\n")
cat("============================================================\n")
cat("ANÁLISIS DE IgG COMPLETADO\n")
cat("============================================================\n\n")

cat("Archivo Excel generado:\n")
cat(output_excel, "\n\n")

cat("Figuras guardadas en:\n")
cat(figures_dir, "\n\n")

cat("Modelos guardados en:\n")
cat(models_dir, "\n\n")

cat("Hojas incluidas en el Excel:\n")

print(
  names(wb)
)

cat("\n============================================================\n")

system("git add R/Curso_temporal_Igg.R data/processed/Curso_temporal_IgG.xlsx results/")
system("git status")
system('git commit -m "Update IgG temporal analysis and results"')
system("git push")

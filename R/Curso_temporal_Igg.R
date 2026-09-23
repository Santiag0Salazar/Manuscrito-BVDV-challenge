# ============================================================
# MANUSCRITO-BVDV-CHALLENGE
# ANÁLISIS TEMPORAL DE IgG ANTI-E2
#
# Curso temporal pre-challenge
# BVDV recombinant E2 vaccine
#
# Figuras:
#   1. Cinética temporal de IgG
#   2. Métricas individuales:
#      Peak, ΔD28-D21, AUC y persistence
#   3. Curvas NLS individuales por grupo
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
  "readr",
  "readxl",
  "openxlsx",
  "nlme",
  "emmeans",
  "broom",
  "patchwork",
  "cowplot",
  "purrr",
  "pracma",
  "scales",
  "splines"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
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
  library(nlme)
  library(emmeans)
  library(broom)
  library(patchwork)
  library(cowplot)
  library(purrr)
  library(pracma)
  library(scales)
  library(splines)
  
})


# ============================================================
# 2. SEMILLA
# ============================================================

set.seed(20260811)


# ============================================================
# 3. LOCALIZAR PROYECTO
# ============================================================

find_project_root <- function() {
  
  current_dir <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
  
  while (
    current_dir != dirname(current_dir)
  ) {
    
    rproj_files <- list.files(
      current_dir,
      pattern = "\\.Rproj$",
      full.names = TRUE
    )
    
    if (length(rproj_files) > 0) {
      return(current_dir)
    }
    
    current_dir <- dirname(current_dir)
  }
  
  stop(
    "No se encontró el archivo .Rproj del proyecto."
  )
}

project_dir <- find_project_root()


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

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  models_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


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
      "No se encontró el archivo:\n",
      input_file
    )
  )
}

cat("\nArchivo de entrada:\n")
cat(input_file, "\n")


# ============================================================
# 6. LECTURA DE DATOS
# ============================================================

datos <- readxl::read_excel(
  input_file,
  sheet = 1
)


# ============================================================
# 7. VALIDACIÓN DE COLUMNAS
# ============================================================

required_columns <- c(
  "Animal",
  "Grupo",
  "Tiempo",
  "Absorbancia"
)

missing_columns <- setdiff(
  required_columns,
  names(datos)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "Faltan las siguientes columnas:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 8. PREPARACIÓN DE DATOS
# ============================================================

datos <- datos %>%
  mutate(
    
    Animal = factor(Animal),
    
    Grupo = case_when(
      
      Grupo %in% c(
        "Comercial",
        "Comm.",
        "Commercial vaccine"
      ) ~ "Commercial",
      
      TRUE ~ as.character(Grupo)
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
  
  mutate(
    
    TiempoF = factor(
      Tiempo,
      levels = c(
        0,
        7,
        14,
        21,
        28,
        35,
        42,
        49,
        56
      )
    ),
    
    TimeIndex = match(
      Tiempo,
      c(
        0,
        7,
        14,
        21,
        28,
        35,
        42,
        49,
        56
      )
    )
    
  )


# ============================================================
# 9. COMPROBACIÓN BÁSICA
# ============================================================

cat("\n================ DATOS =================\n")

print(
  datos %>%
    count(Grupo, Tiempo)
)

cat("\nNúmero de animales:\n")

print(
  datos %>%
    distinct(Animal, Grupo) %>%
    count(Grupo)
)


# ============================================================
# 10. DATOS VACUNADOS
# ============================================================

datos_vac <- datos %>%
  filter(
    Grupo != "Control"
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


# Etiquetas utilizadas en los gráficos
group_labels <- c(
  
  "Control" = "Control",
  
  "25ug" = "25 µg",
  
  "50ug" = "50 µg",
  
  "100ug" = "100 µg",
  
  "Commercial" = "Comm."
  
)


# ============================================================
# 12. TEMA
# ============================================================

theme_nature <- function(base_size = 10) {
  
  cowplot::theme_cowplot(
    font_size = base_size,
    font_family = "sans"
  ) +
    
    theme(
      
      plot.title = element_text(
        face = "plain",
        size = base_size,
        hjust = 0.5
      ),
      
      axis.title = element_text(
        size = base_size
      ),
      
      axis.text = element_text(
        size = base_size,
        color = "grey15"
      ),
      
      legend.position = "none",
      
      panel.grid.major.y = element_line(
        color = "grey90",
        linewidth = 0.25
      ),
      
      panel.grid.major.x = element_blank()
      
    )
}


# ============================================================
# 13. FUNCIÓN PARA GUARDAR FIGURAS
# ============================================================

save_nature <- function(
    plot,
    filename,
    width = 7,
    height = 4
) {
  
  ggsave(
    filename = file.path(
      figures_dir,
      paste0(filename, ".pdf")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )
  
  ggsave(
    filename = file.path(
      figures_dir,
      paste0(filename, ".tiff")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    filename = file.path(
      figures_dir,
      paste0(filename, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
}


# ============================================================
# 14. MEDIA Y 95% CI
# ============================================================

mean_ci <- datos %>%
  
  group_by(
    Grupo,
    Tiempo
  ) %>%
  
  summarise(
    
    n = sum(
      !is.na(Absorbancia)
    ),
    
    mean = mean(
      Absorbancia,
      na.rm = TRUE
    ),
    
    sd = sd(
      Absorbancia,
      na.rm = TRUE
    ),
    
    se = sd / sqrt(n),
    
    lower = mean - qt(
      0.975,
      df = n - 1
    ) * se,
    
    upper = mean + qt(
      0.975,
      df = n - 1
    ) * se,
    
    .groups = "drop"
    
  )


# ============================================================
# 15. MÉTRICAS INDIVIDUALES POR ANIMAL
# ============================================================

animal_metrics <- datos %>%
  
  group_by(
    Animal,
    Grupo
  ) %>%
  
  arrange(
    Tiempo,
    .by_group = TRUE
  ) %>%
  
  summarise(
    
    peak = max(
      Absorbancia,
      na.rm = TRUE
    ),
    
    D21 = ifelse(
      any(Tiempo == 21),
      Absorbancia[Tiempo == 21][1],
      NA_real_
    ),
    
    D28 = ifelse(
      any(Tiempo == 28),
      Absorbancia[Tiempo == 28][1],
      NA_real_
    ),
    
    D56 = ifelse(
      any(Tiempo == 56),
      Absorbancia[Tiempo == 56][1],
      NA_real_
    ),
    
    AUC = if (
      sum(!is.na(Absorbancia)) >= 2
    ) {
      
      pracma::trapz(
        Tiempo[!is.na(Absorbancia)],
        Absorbancia[!is.na(Absorbancia)]
      )
      
    } else {
      NA_real_
    },
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    
    delta_D28_D21 =
      D28 - D21,
    
    persistence =
      D56 / peak
    
  )


# ============================================================
# 16. RESUMEN DE MÉTRICAS
# ============================================================

metric_summary <- animal_metrics %>%
  
  group_by(
    Grupo
  ) %>%
  
  summarise(
    
    n = n(),
    
    peak_mean = mean(
      peak,
      na.rm = TRUE
    ),
    
    peak_sd = sd(
      peak,
      na.rm = TRUE
    ),
    
    delta_mean = mean(
      delta_D28_D21,
      na.rm = TRUE
    ),
    
    delta_sd = sd(
      delta_D28_D21,
      na.rm = TRUE
    ),
    
    AUC_mean = mean(
      AUC,
      na.rm = TRUE
    ),
    
    AUC_sd = sd(
      AUC,
      na.rm = TRUE
    ),
    
    persistence_mean = mean(
      persistence,
      na.rm = TRUE
    ),
    
    persistence_sd = sd(
      persistence,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  )


# ============================================================
# 17. FUNCIÓN MODELO MIXTO
# ============================================================

fit_lme <- function(
    fixed_formula,
    data,
    method = "REML",
    correlation = NULL
) {
  
  tryCatch(
    
    {
      
      if (is.null(correlation)) {
        
        nlme::lme(
          
          fixed = fixed_formula,
          
          random = ~1 | Animal,
          
          data = data,
          
          method = method,
          
          na.action = na.omit,
          
          keep.data = TRUE
          
        )
        
      } else {
        
        nlme::lme(
          
          fixed = fixed_formula,
          
          random = ~1 | Animal,
          
          data = data,
          
          method = method,
          
          correlation = correlation,
          
          na.action = na.omit,
          
          keep.data = TRUE
          
        )
        
      }
      
    },
    
    error = function(e) {
      
      message(
        "Modelo no pudo ajustarse: ",
        e$message
      )
      
      NULL
      
    }
    
  )
}


# ============================================================
# 18. MODELOS MIXTOS CANDIDATOS
# ============================================================

model_null <- fit_lme(
  Absorbancia ~ 1,
  datos
)

model_group <- fit_lme(
  Absorbancia ~ Grupo,
  datos
)

model_time <- fit_lme(
  Absorbancia ~ TiempoF,
  datos
)

model_group_time <- fit_lme(
  Absorbancia ~ Grupo + TiempoF,
  datos
)

model_interaction <- fit_lme(
  Absorbancia ~ Grupo * TiempoF,
  datos
)


# ============================================================
# 19. COMPARACIÓN DE MODELOS
# ============================================================

model_list <- list(
  
  Null = model_null,
  
  Group = model_group,
  
  Time = model_time,
  
  Group_Time = model_group_time,
  
  Group_Time_Interaction = model_interaction
  
)

model_list <- model_list[
  !vapply(
    model_list,
    is.null,
    logical(1)
  )
]


model_comparison <- bind_rows(
  
  lapply(
    names(model_list),
    function(x) {
      
      m <- model_list[[x]]
      
      data.frame(
        
        Model = x,
        
        AIC = AIC(m),
        
        BIC = BIC(m),
        
        logLik = as.numeric(
          logLik(m)
        )
        
      )
      
    }
  )
  
) %>%
  
  arrange(AIC)


cat(
  "\n================ COMPARACIÓN DE MODELOS ================\n"
)

print(model_comparison)


# ============================================================
# 20. MODELO FINAL
# ============================================================

final_time_model <- model_interaction


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
# 22. EMMEANS
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
# 23. COMPARACIONES POST-HOC ENTRE GRUPOS
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
# 24. ANÁLISIS DE MÉTRICAS
# ============================================================

metric_tests <- list()


metrics_to_test <- c(
  
  "peak",
  
  "delta_D28_D21",
  
  "AUC",
  
  "persistence"
  
)


for (
  metric_name in metrics_to_test
) {
  
  df_metric <- animal_metrics %>%
    
    select(
      Grupo,
      all_of(metric_name)
    ) %>%
    
    rename(
      Value = all_of(metric_name)
    ) %>%
    
    filter(
      !is.na(Value)
    )
  
  if (
    nrow(df_metric) > 0
  ) {
    
    fit <- tryCatch(
      
      aov(
        Value ~ Grupo,
        data = df_metric
      ),
      
      error = function(e) NULL
      
    )
    
    if (!is.null(fit)) {
      
      metric_tests[[metric_name]] <- list(
        
        anova = broom::tidy(fit),
        
        tukey = broom::tidy(
          TukeyHSD(
            fit,
            "Grupo"
          )
        )
        
      )
      
    }
    
  }
  
}


# ============================================================
# 25. FUNCIÓN PARA BOXPLOTS
# ESTILO AUC SOLICITADO
# ============================================================

generate_metric_boxplot <- function(
    
  df,
  
  metric,
  
  y_label,
  
  title = NULL,
  
  y_breaks = waiver()
  
) {
  
  df_plot <- df %>%
    
    select(
      Grupo,
      all_of(metric)
    ) %>%
    
    rename(
      Value = all_of(metric)
    ) %>%
    
    filter(
      !is.na(Value)
    ) %>%
    
    mutate(
      
      Grupo = factor(
        Grupo,
        levels = c(
          "Control",
          "25ug",
          "50ug",
          "100ug",
          "Commercial"
        )
      )
      
    )
  
  
  # ----------------------------------------------------------
  # ANOVA
  # ----------------------------------------------------------
  
  fit_aov <- aov(
    Value ~ Grupo,
    data = df_plot
  )
  
  
  # ----------------------------------------------------------
  # TUKEY
  # ----------------------------------------------------------
  
  tukey_out <- TukeyHSD(
    fit_aov,
    "Grupo"
  )$Grupo
  
  
  tukey_df <- data.frame(
    
    comparisons = rownames(
      tukey_out
    ),
    
    p.value = tukey_out[
      ,
      "p adj"
    ],
    
    row.names = NULL
    
  )
  
  
  # ----------------------------------------------------------
  # RANGO Y
  # ----------------------------------------------------------
  
  max_y <- max(
    df_plot$Value,
    na.rm = TRUE
  )
  
  min_y <- min(
    df_plot$Value,
    na.rm = TRUE
  )
  
  range_y <- max_y - min_y
  
  if (
    range_y == 0 ||
    !is.finite(range_y)
  ) {
    
    range_y <- 1
    
  }
  
  
  # ----------------------------------------------------------
  # COMPARACIONES
  #
  # 25ug vs Commercial
  # 50ug vs Commercial
  # 100ug vs Commercial
  # ----------------------------------------------------------
  
  comp_metadata <- list(
    
    list(
      x1 = 2,
      x2 = 5,
      row = 1,
      pair_pattern =
        "Commercial-25ug|25ug-Commercial"
    ),
    
    list(
      x1 = 3,
      x2 = 5,
      row = 2,
      pair_pattern =
        "Commercial-50ug|50ug-Commercial"
    ),
    
    list(
      x1 = 4,
      x2 = 5,
      row = 3,
      pair_pattern =
        "Commercial-100ug|100ug-Commercial"
    )
    
  )
  
  
  lines_and_labels <- purrr::map_df(
    
    comp_metadata,
    
    function(m) {
      
      p_val <- tukey_df %>%
        
        filter(
          grepl(
            m$pair_pattern,
            comparisons
          )
        ) %>%
        
        pull(
          p.value
        )
      
      
      if (
        length(p_val) == 0
      ) {
        
        p_val <- 1
        
      }
      
      
      p_val <- p_val[1]
      
      
      p_formatted <- if (
        p_val < 0.001
      ) {
        
        "p < 0.001"
        
      } else if (
        p_val >= 0.05
      ) {
        
        "ns"
        
      } else {
        
        paste0(
          "p = ",
          sprintf(
            "%.3f",
            p_val
          )
        )
        
      }
      
      
      y_bar <- max_y +
        range_y *
        0.10 *
        m$row
      
      
      data.frame(
        
        x = m$x1,
        
        xend = m$x2,
        
        y = y_bar,
        
        x_text =
          (m$x1 + m$x2) / 2,
        
        y_text =
          y_bar +
          range_y * 0.025,
        
        p_label =
          p_formatted
        
      )
      
    }
    
  )
  
  
  ylim_max <- max(
    lines_and_labels$y_text
  ) +
    range_y * 0.08
  
  
  # ----------------------------------------------------------
  # GRÁFICO
  # ----------------------------------------------------------
  
  p <- ggplot(
    
    df_plot,
    
    aes(
      x = Grupo,
      y = Value,
      color = Grupo
    )
    
  ) +
    
    geom_boxplot(
      
      aes(
        fill = Grupo
      ),
      
      width = 0.45,
      
      outlier.shape = NA,
      
      alpha = 0.15,
      
      linewidth = 0.4
      
    ) +
    
    geom_jitter(
      
      width = 0.10,
      
      height = 0,
      
      size = 1.4,
      
      alpha = 0.85
      
    ) +
    
    
    # Barra horizontal
    geom_segment(
      
      data = lines_and_labels,
      
      aes(
        x = x,
        xend = xend,
        y = y,
        yend = y
      ),
      
      color = "grey35",
      
      linewidth = 0.35,
      
      inherit.aes = FALSE
      
    ) +
    
    
    # Barra vertical izquierda
    geom_segment(
      
      data = lines_and_labels,
      
      aes(
        x = x,
        xend = x,
        y = y,
        yend =
          y -
          range_y * 0.015
      ),
      
      color = "grey35",
      
      linewidth = 0.35,
      
      inherit.aes = FALSE
      
    ) +
    
    
    # Barra vertical derecha
    geom_segment(
      
      data = lines_and_labels,
      
      aes(
        x = xend,
        xend = xend,
        y = y,
        yend =
          y -
          range_y * 0.015
      ),
      
      color = "grey35",
      
      linewidth = 0.35,
      
      inherit.aes = FALSE
      
    ) +
    
    
    # P-values
    geom_text(
      
      data = lines_and_labels,
      
      aes(
        x = x_text,
        y = y_text,
        label = p_label
      ),
      
      color = "black",
      
      size = 2.2,
      
      vjust = 0,
      
      fontface = "plain",
      
      inherit.aes = FALSE
      
    ) +
    
    
    scale_color_manual(
      
      values = palette_nature,
      
      labels = group_labels,
      
      drop = FALSE
      
    ) +
    
    
    scale_fill_manual(
      
      values = palette_nature,
      
      labels = group_labels,
      
      drop = FALSE
      
    ) +
    
    
    coord_cartesian(
      
      ylim = c(
        
        min_y -
          range_y * 0.05,
        
        ylim_max
        
      )
      
    ) +
    
    
    labs(
      
      title = title,
      
      x = NULL,
      
      y = y_label
      
    ) +
    
    
    theme_nature(
      11
    ) +
    
    
    theme(
      
      legend.position = "bottom",
      
      legend.justification = "center"
      
    )
  
  
  # ----------------------------------------------------------
  # ESCALA Y
  # ----------------------------------------------------------
  
  if (
    !is.null(y_breaks)
  ) {
    
    p <- p +
      
      scale_y_continuous(
        breaks = y_breaks
      )
    
  }
  
  
  return(p)
  
}


# ============================================================
# 26. FIGURA 2A — PEAK
# ============================================================

fig_peak <- generate_metric_boxplot(
  
  df = animal_metrics,
  
  metric = "peak",
  
  y_label =
    "Peak anti-E2 IgG absorbance",
  
  title =
    "Peak anti-E2 IgG"
  
)


# ============================================================
# 27. FIGURA 2B — ΔD28-D21
# ============================================================

fig_delta <- generate_metric_boxplot(
  
  df = animal_metrics,
  
  metric = "delta_D28_D21",
  
  y_label =
    expression(
      Delta * "D28–D21 anti-E2 IgG"
    ),
  
  title =
    "Net increase after booster"
  
)


# ============================================================
# 28. FIGURA 2C — AUC
# ============================================================

fig_auc <- generate_metric_boxplot(
  
  df = animal_metrics,
  
  metric = "AUC",
  
  y_label =
    expression(
      "Area Under the Curve"
    ),
  
  title =
    "Area Under the Titer Curve"
  
)


# ============================================================
# 29. FIGURA 2D — PERSISTENCE
# ============================================================

fig_persistence <- generate_metric_boxplot(
  
  df = animal_metrics,
  
  metric = "persistence",
  
  y_label =
    "Persistence (D56 / peak)",
  
  title =
    "Antibody persistence"
  
)


# ============================================================
# 30. FIGURA 2 COMPLETA
# ============================================================

fig2_final <- (
  
  fig_peak +
    
    fig_delta +
    
    fig_auc +
    
    fig_persistence
  
) +
  
  plot_layout(
    ncol = 4
  )


# ============================================================
# 31. GUARDAR FIGURA 2
# ============================================================

save_nature(
  
  fig2_final,
  
  "IgG_Figure_2_individual_metrics",
  
  width = 12,
  
  height = 3.8
  
)


# ============================================================
# 32. FIGURA 1
# CINÉTICA TEMPORAL
# ============================================================

fig1_final <- ggplot(
  
  mean_ci,
  
  aes(
    x = Tiempo,
    y = mean,
    color = Grupo,
    fill = Grupo
  )
  
) +
  
  geom_ribbon(
    
    aes(
      ymin = lower,
      ymax = upper
    ),
    
    alpha = 0.12,
    
    color = NA
    
  ) +
  
  geom_line(
    
    linewidth = 0.8
    
  ) +
  
  geom_point(
    
    size = 2,
    
    alpha = 0.85
    
  ) +
  
  geom_vline(
    
    xintercept = 21,
    
    linetype = "dashed",
    
    linewidth = 0.35,
    
    color = "grey40",
    
    alpha = 0.5
    
  ) +
  
  annotate(
    
    "text",
    
    x = 22,
    
    y = max(
      mean_ci$upper,
      na.rm = TRUE
    ),
    
    label = "Booster",
    
    hjust = 0,
    
    vjust = 1,
    
    size = 3.5,
    
    fontface = "italic",
    
    color = "grey30"
    
  ) +
  
  scale_color_manual(
    
    values = palette_nature,
    
    labels = group_labels,
    
    drop = FALSE
    
  ) +
  
  scale_fill_manual(
    
    values = palette_nature,
    
    labels = group_labels,
    
    drop = FALSE
    
  ) +
  
  scale_x_continuous(
    
    breaks = c(
      0,
      7,
      14,
      21,
      28,
      35,
      42,
      49,
      56
    ),
    
    expand =
      expansion(
        mult = c(
          0.02,
          0.04
        )
      )
    
  ) +
  
  labs(
    
    x =
      "Time after first immunization (days)",
    
    y =
      "Anti-E2 IgG (Absorbance, 450 nm)",
    
    color =
      "Vaccine Group",
    
    fill =
      "Vaccine Group"
    
  ) +
  
  theme_nature(
    10
  ) +
  
  theme(
    
    legend.position = "top",
    
    legend.justification = "center"
    
  )


# ============================================================
# 33. GUARDAR FIGURA 1
# ============================================================

save_nature(
  
  fig1_final,
  
  "IgG_Figure_1_temporal_kinetics",
  
  width = 7.2,
  
  height = 3.2
  
)


# ============================================================
# 34. NLS POR GRUPO
#
# Logistic
# Gompertz
# Weibull
#
# Control NO se ajusta mediante NLS.
# Se representa mediante la trayectoria empírica.
# ============================================================

group_means <- datos_vac %>%
  
  group_by(
    Grupo,
    Tiempo
  ) %>%
  
  summarise(
    
    Absorbancia =
      mean(
        Absorbancia,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )


fit_nls_models <- function(df) {
  
  y_min <- min(
    df$Absorbancia,
    na.rm = TRUE
  )
  
  y_max <- max(
    df$Absorbancia,
    na.rm = TRUE
  )
  
  amp <- max(
    y_max - y_min,
    0.1
  )
  
  
  safe_nls <- function(
    
    formula,
    
    start,
    
    lower,
    
    upper,
    
    model_name
    
  ) {
    
    fit <- tryCatch(
      
      suppressWarnings(
        
        nls(
          
          formula,
          
          data = df,
          
          start = start,
          
          algorithm = "port",
          
          lower = lower,
          
          upper = upper,
          
          control =
            nls.control(
              maxiter = 500,
              warnOnly = TRUE
            )
          
        )
        
      ),
      
      error = function(e) {
        
        NULL
        
      }
      
    )
    
    
    if (
      is.null(fit)
    ) {
      
      return(NULL)
      
    }
    
    
    list(
      
      model =
        model_name,
      
      fit =
        fit,
      
      AIC =
        AIC(fit)
      
    )
    
  }
  
  
  list(
    
    Logistic = safe_nls(
      
      Absorbancia ~
        
        Bottom +
        
        (Top - Bottom) /
        
        (
          1 +
            exp(
              -(Tiempo - T50) / s
            )
        ),
      
      start = list(
        
        Bottom = y_min,
        
        Top = y_max,
        
        T50 = 18,
        
        s = 3
        
      ),
      
      lower = c(
        
        Bottom = 0,
        
        Top = 0.3,
        
        T50 = 0,
        
        s = 0.2
        
      ),
      
      upper = c(
        
        Bottom = 1.0,
        
        Top = 3.0,
        
        T50 = 56,
        
        s = 20
        
      ),
      
      model_name = "Logistic"
      
    ),
    
    
    Gompertz = safe_nls(
      
      Absorbancia ~
        
        Bottom +
        
        A *
        
        exp(
          -exp(
            -k * (Tiempo - Ti)
          )
        ),
      
      start = list(
        
        Bottom = y_min,
        
        A = amp,
        
        k = 0.25,
        
        Ti = 14
        
      ),
      
      lower = c(
        
        Bottom = 0,
        
        A = 0.1,
        
        k = 0.001,
        
        Ti = 0
        
      ),
      
      upper = c(
        
        Bottom = 1.0,
        
        A = 3.0,
        
        k = 2.0,
        
        Ti = 56
        
      ),
      
      model_name = "Gompertz"
      
    ),
    
    
    Weibull = safe_nls(
      
      Absorbancia ~
        
        Bottom +
        
        A *
        
        (
          1 -
            exp(
              -(Tiempo / lambda)^beta
            )
        ),
      
      start = list(
        
        Bottom = y_min,
        
        A = amp,
        
        lambda = 16,
        
        beta = 3
        
      ),
      
      lower = c(
        
        Bottom = 0,
        
        A = 0.1,
        
        lambda = 1,
        
        beta = 0.2
        
      ),
      
      upper = c(
        
        Bottom = 1.0,
        
        A = 3.0,
        
        lambda = 80,
        
        beta = 20
        
      ),
      
      model_name = "Weibull"
      
    )
    
  )
  
}


nls_fits <- lapply(
  
  split(
    group_means,
    group_means$Grupo
  ),
  
  fit_nls_models
  
)


# ============================================================
# 35. SELECCIONAR MEJOR MODELO POR GRUPO
# ============================================================

nls_summary <- bind_rows(
  
  lapply(
    
    names(nls_fits),
    
    function(grp) {
      
      bind_rows(
        
        lapply(
          
          nls_fits[[grp]],
          
          function(x) {
            
            if (
              is.null(x)
            ) {
              
              return(NULL)
              
            }
            
            data.frame(
              
              Grupo = grp,
              
              Modelo = x$model,
              
              AIC = x$AIC
              
            )
            
          }
          
        )
        
      )
      
    }
    
  )
  
) %>%
  
  group_by(
    Grupo
  ) %>%
  
  slice_min(
    
    order_by = AIC,
    
    n = 1,
    
    with_ties = FALSE
    
  ) %>%
  
  ungroup()


cat(
  "\n================ MEJOR MODELO NLS POR GRUPO ================\n"
)

print(
  nls_summary
)


# ============================================================
# 36. PREDICCIONES NLS
# ============================================================

x_pred <- seq(
  
  min(
    datos$Tiempo,
    na.rm = TRUE
  ),
  
  max(
    datos$Tiempo,
    na.rm = TRUE
  ),
  
  length.out = 300
  
)


nls_predictions <- bind_rows(
  
  lapply(
    
    seq_len(
      nrow(nls_summary)
    ),
    
    function(i) {
      
      grp <-
        as.character(
          nls_summary$Grupo[i]
        )
      
      model_name <-
        as.character(
          nls_summary$Modelo[i]
        )
      
      fit <-
        nls_fits[[grp]][[model_name]]$fit
      
      
      nd <- data.frame(
        Tiempo = x_pred
      )
      
      
      data.frame(
        
        Grupo = grp,
        
        Modelo = model_name,
        
        Tiempo = x_pred,
        
        Absorbancia =
          as.numeric(
            predict(
              fit,
              nd
            )
          )
        
      )
      
    }
    
  )
  
)


# ============================================================
# 37. FIGURA 3
# CURVAS NLS INDIVIDUALES POR GRUPO
# ============================================================

group_levels <- c(
  
  "Control",
  
  "25ug",
  
  "50ug",
  
  "100ug",
  
  "Commercial"
  
)


datos$Grupo <-
  factor(
    datos$Grupo,
    levels = group_levels
  )


mean_ci$Grupo <-
  factor(
    mean_ci$Grupo,
    levels = group_levels
  )


nls_predictions$Grupo <-
  factor(
    nls_predictions$Grupo,
    levels = group_levels
  )


fig3_final_con_control <- ggplot() +
  
  # ----------------------------------------------------------
# Datos individuales
# ----------------------------------------------------------

geom_jitter(
  
  data = datos,
  
  aes(
    x = Tiempo,
    y = Absorbancia,
    color = Grupo
  ),
  
  size = 0.9,
  
  alpha = 0.15,
  
  width = 0.6,
  
  height = 0
  
) +
  
  
  # ----------------------------------------------------------
# Control: trayectoria empírica
# ----------------------------------------------------------

geom_line(
  
  data =
    mean_ci %>%
    filter(
      Grupo == "Control"
    ),
  
  aes(
    x = Tiempo,
    y = mean,
    color = Grupo
  ),
  
  linewidth = 0.6,
  
  linetype = "twodash",
  
  alpha = 0.7
  
) +
  
  
  # ----------------------------------------------------------
# Curvas NLS
# ----------------------------------------------------------

geom_line(
  
  data = nls_predictions,
  
  aes(
    x = Tiempo,
    y = Absorbancia,
    color = Grupo
  ),
  
  linewidth = 0.95,
  
  alpha = 0.95
  
) +
  
  
  # ----------------------------------------------------------
# Medias
# ----------------------------------------------------------

geom_point(
  
  data = mean_ci,
  
  aes(
    x = Tiempo,
    y = mean,
    color = Grupo
  ),
  
  size = 2.0,
  
  alpha = 0.75,
  
  shape = 16
  
) +
  
  
  # ----------------------------------------------------------
# Booster
# ----------------------------------------------------------

geom_vline(
  
  xintercept = 21,
  
  linetype = "dashed",
  
  linewidth = 0.35,
  
  color = "grey40",
  
  alpha = 0.5
  
) +
  
  
  annotate(
    
    "text",
    
    x = 22,
    
    y = 2.0,
    
    label = "Booster",
    
    hjust = 0,
    
    vjust = 1,
    
    size = 4,
    
    fontface = "italic",
    
    color = "grey30"
    
  ) +
  
  
  scale_color_manual(
    
    values = palette_nature,
    
    labels = group_labels,
    
    drop = FALSE
    
  ) +
  
  
  scale_x_continuous(
    
    breaks = c(
      
      0,
      
      7,
      
      14,
      
      21,
      
      28,
      
      35,
      
      42,
      
      49,
      
      56
      
    ),
    
    expand =
      expansion(
        mult = c(
          0.02,
          0.04
        )
      )
    
  ) +
  
  
  scale_y_continuous(
    
    breaks =
      seq(
        0,
        2.0,
        0.5
      ),
    
    expand =
      expansion(
        mult = c(
          0,
          0.03
        )
      )
    
  ) +
  
  
  coord_cartesian(
    
    ylim = c(
      0,
      2.15
    )
    
  ) +
  
  
  labs(
    
    x =
      "Time after first immunization (days)",
    
    y =
      "Anti-E2 IgG (Absorbance, 450 nm)",
    
    color =
      "Vaccine Group"
    
  ) +
  
  
  theme_nature(
    10
  ) +
  
  
  theme(
    
    legend.position = "top",
    
    legend.justification = "center"
    
  )


# ============================================================
# 38. GUARDAR FIGURA 3
# ============================================================

save_nature(
  
  fig3_final_con_control,
  
  "Curso temporal_anti_E2_Igg",
  
  width = 7.2,
  
  height = 3.2
  
)


# ============================================================
# 39. PARÁMETROS NLS
# ============================================================

nls_parameters_final <- bind_rows(
  
  lapply(
    
    seq_len(
      nrow(nls_summary)
    ),
    
    function(i) {
      
      grp <-
        as.character(
          nls_summary$Grupo[i]
        )
      
      model_name <-
        as.character(
          nls_summary$Modelo[i]
        )
      
      fit <-
        nls_fits[[grp]][[model_name]]$fit
      
      
      broom::tidy(
        fit
      ) %>%
        
        mutate(
          
          Grupo = grp,
          
          Modelo = model_name,
          
          AIC = AIC(fit),
          
          .before = 1
          
        )
      
    }
    
  )
  
)


# ============================================================
# 40. EXPORTAR INFORMACIÓN DE DATOS
# ============================================================

analysis_info <- data.frame(
  
  Item = c(
    
    "Analysis",
    
    "Input file",
    
    "Number of animals",
    
    "Number of observations",
    
    "Groups",
    
    "Time points",
    
    "Longitudinal model",
    
    "Random effect",
    
    "Multiple comparison",
    
    "NLS models"
    
  ),
  
  Value = c(
    
    "Temporal anti-E2 IgG analysis",
    
    basename(input_file),
    
    n_distinct(
      datos$Animal
    ),
    
    nrow(datos),
    
    paste(
      levels(datos$Grupo),
      collapse = ", "
    ),
    
    paste(
      sort(
        unique(
          datos$Tiempo
        )
      ),
      collapse = ", "
    ),
    
    "Linear mixed-effects model",
    
    "Animal",
    
    "Tukey",
    
    "Logistic, Gompertz, Weibull"
    
  )
  
)


# ============================================================
# 41. GUARDAR INFO DEL ANÁLISIS
# ============================================================

write.csv(
  
  analysis_info,
  
  file.path(
    tables_dir,
    "IgG_analysis_info.csv"
  ),
  
  row.names = FALSE
  
)


# ============================================================
# 42. EXPORTAR MODELOS
# ============================================================

if (
  !is.null(final_time_model)
) {
  
  saveRDS(
    
    final_time_model,
    
    file.path(
      models_dir,
      "IgG_modelo_mixto_tiempo_REML.rds"
    )
    
  )
  
}


# ============================================================
# 43. GUARDAR MODELOS NLS
# ============================================================

if (
  !is.null(nls_fits)
) {
  
  best_nls_models <- lapply(
    
    names(nls_fits),
    
    function(grp) {
      
      model_name <-
        nls_summary %>%
        filter(
          Grupo == grp
        ) %>%
        pull(
          Modelo
        )
      
      if (
        length(model_name) == 0
      ) {
        
        return(NULL)
        
      }
      
      nls_fits[[grp]][[
        model_name[1]
      ]]
      
    }
    
  )
  
  names(
    best_nls_models
  ) <-
    names(nls_fits)
  
  
  saveRDS(
    
    best_nls_models,
    
    file.path(
      models_dir,
      "IgG_NLS_modelos_por_grupo.rds"
    )
    
  )
  
}


# ============================================================
# 44. EXPORTAR RESULTADOS A UN SOLO EXCEL
# ============================================================

output_excel <- file.path(
  
  tables_dir,
  
  "Curso_temporal_IgG_resultados.xlsx"
  
)


wb <- createWorkbook()


add_sheet <- function(
    wb,
    sheet_name,
    data
) {
  
  addWorksheet(
    wb,
    sheet_name
  )
  
  if (
    is.null(data) ||
    nrow(data) == 0
  ) {
    
    writeData(
      wb,
      sheet = sheet_name,
      x = data.frame(
        Information =
          "No data available"
      )
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
  
}


# ------------------------------------------------------------
# README
# ------------------------------------------------------------

readme_table <- data.frame(
  
  Section = c(
    
    "Input",
    
    "Temporal analysis",
    
    "Metric analysis",
    
    "Figure 1",
    
    "Figure 2",
    
    "Figure 3",
    
    "NLS"
    
  ),
  
  Description = c(
    
    "Original processed IgG dataset",
    
    "Linear mixed-effects model with animal random intercept",
    
    "Peak, ΔD28-D21, AUC and persistence",
    
    "Temporal anti-E2 IgG kinetics",
    
    "Individual metric boxplots with Tukey comparisons",
    
    "Individual NLS curves by vaccine group",
    
    "Logistic, Gompertz and Weibull selected by AIC"
    
  )
  
)


add_sheet(
  wb,
  "README",
  readme_table
)


# ------------------------------------------------------------
# Datos
# ------------------------------------------------------------

add_sheet(
  wb,
  "Datos",
  datos
)


# ------------------------------------------------------------
# Resumen temporal
# ------------------------------------------------------------

add_sheet(
  wb,
  "Resumen",
  mean_ci
)


# ------------------------------------------------------------
# Métricas por animal
# ------------------------------------------------------------

add_sheet(
  wb,
  "Metricas_animal",
  animal_metrics
)


# ------------------------------------------------------------
# Resumen métricas
# ------------------------------------------------------------

add_sheet(
  wb,
  "Metricas_resumen",
  metric_summary
)


# ------------------------------------------------------------
# Modelo mixto
# ------------------------------------------------------------

add_sheet(
  wb,
  "Modelo_mixto",
  model_comparison
)


# ------------------------------------------------------------
# ANOVA modelo final
# ------------------------------------------------------------

add_sheet(
  wb,
  "ANOVA_modelo",
  anova_final
)


# ------------------------------------------------------------
# EMMEANS
# ------------------------------------------------------------

add_sheet(
  wb,
  "EMMmeans",
  emmeans_time
)


# ------------------------------------------------------------
# Comparaciones tiempo
# ------------------------------------------------------------

add_sheet(
  wb,
  "Comparaciones_tiempo",
  pairwise_time
)


# ------------------------------------------------------------
# NLS mejor modelo por grupo
# ------------------------------------------------------------

add_sheet(
  
  wb,
  
  "NLS_mejor_modelo_grupo",
  
  nls_summary
  
)


# ------------------------------------------------------------
# NLS parámetros
# ------------------------------------------------------------

add_sheet(
  
  wb,
  
  "NLS_parametros_grupo",
  
  nls_parameters_final
  
)


# ------------------------------------------------------------
# NLS predicciones
# ------------------------------------------------------------

add_sheet(
  
  wb,
  
  "NLS_predicciones",
  
  nls_predictions
  
)


# ------------------------------------------------------------
# Información análisis
# ------------------------------------------------------------

add_sheet(
  
  wb,
  
  "Analysis_info",
  
  analysis_info
  
)


# ============================================================
# 45. GUARDAR EXCEL
# ============================================================

saveWorkbook(
  
  wb,
  
  output_excel,
  
  overwrite = TRUE
  
)


# ============================================================
# 46. MENSAJE FINAL
# ============================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "ANÁLISIS COMPLETADO CORRECTAMENTE\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nResultados:\n"
)

cat(
  output_excel,
  "\n"
)

cat(
  "\nFiguras guardadas en:\n"
)

cat(
  figures_dir,
  "\n"
)

cat(
  "\nModelos guardados en:\n"
)

cat(
  models_dir,
  "\n"
)

cat(
  "\nFigura 2: boxplots + jitter + Tukey\n"
)

cat(
  "Figura 3: NLS por grupo\n"
)

cat(
  "\n============================================================\n"
)
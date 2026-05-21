# ==============================================================================
# ANÁLISIS DE FRECUENCIAS DE CAUDALES MÁXIMOS — CFE Colombia
# ==============================================================================
# Estación   : Río Guachicono (Cauca, Colombia)
# Variable   : Caudal máximo mensual multianual (m³/s)
# Archivo    : D:/R/Version_Final_CFE_Colombia/Qmax.xlsx
# Estructura : Año | Enero | Febrero | ... | Diciembre  (13 columnas)
#
# Autor      : Mauricio Javier Victoria Niño
# Filiación  : Investigador independiente · Cali, Colombia
# Contacto   : hidratecsa@gmail.com
# ORCID      : 0009-0003-4328-5691
#
# Versión    : 1.0.0  (2026-05-20)
# Licencia   : CC BY 4.0  https://creativecommons.org/licenses/by/4.0/
#
# CÓMO USAR :
#   1. Colocar Qmax.xlsx en D:/R/Version_Final_CFE_Colombia/
#   2. La carpeta outputs/ se crea automáticamente.
#   3. Ejecutar el script completo (Ctrl+Alt+R en RStudio).
#
# MÓDULOS:
#   1. Detección y caracterización de no estacionariedad
#      1.1 Análisis exploratorio visual (panel 4-en-1)
#      1.2 Test Ljung-Box + batería de 6 tests de estacionariedad
#      1.3 Tendencia por época hidrológica
#      1.4 Frecuencia de extremos en ventanas móviles
#      1.5 Análisis wavelet multi-escala (Morlet, n.sim = 1000)
#      1.6 Clasificación final con criterio documentado
#
#   2. Análisis de frecuencias de caudales máximos
#      2.1 Estacionario: LP3, GEV, GUM, LNO, PE3, NOR
#          · Selección por AICc (Burnham & Anderson 2002)
#          · Posición de trazado Cunnane
#          · IC 90 % bootstrap paramétrico GEV (B = 1000)
#      2.2 No estacionario: GAMLSS (6 modelos)
#          · Validación Cox-Snell: fitted() + KS + Filliben
#          · IC 95 % bootstrap de residuos (B = 500)
#      2.3 Informe final con sessionInfo()
#
# REFERENCIAS:
#   Burnham & Anderson (2002) Model Selection and Multimodel Inference. Springer.
#   Rigby & Stasinopoulos (2005) Appl. Stat. 54(3):507-554.
#   Villarini et al. (2009) Adv. Water Resour. 32:1255-1266.
#   Yue & Wang (2002) Adv. Water Resour. 25(3):325-333.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURACIÓN
# ==============================================================================

# 0.1 Paquetes -----------------------------------------------------------------
required_packages <- c(
  "readxl", "dplyr", "lubridate", "ggplot2", "tidyr", "gridExtra",
  "gamlss", "gamlss.dist", "gamlss.add",
  "ismev", "extRemes", "fitdistrplus", "MASS", "scales",
  "trend", "Kendall", "moments", "evd", "WaveletComp",
  "zoo", "tseries", "lmtest", "janitor"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}
library(grid)

# 0.2 Funciones auxiliares -----------------------------------------------------

encontrar_idx <- function(valor, vector, tol = 0.001) {
  idx <- which(abs(vector - valor) < tol)
  if (length(idx) > 0) idx[1] else NA
}

calcular_u <- function(x) {
  n <- length(x); u <- rep(NA, n)
  for (k in 2:n) {
    ti <- 0
    for (j in 1:k) ti <- ti + sum(x[j] > x[1:j])
    u[k] <- (ti - k*(k-1)/4) / sqrt(k*(k-1)*(2*k+5)/72)
  }
  u
}

posicion_trazado <- function(n, metodo = "cunnane") {
  i <- seq_len(n)
  switch(metodo,
    weibull    = i / (n + 1),
    cunnane    = (i - 0.40) / (n + 0.20),
    gringorten = (i - 0.44) / (n + 0.12),
    hazen      = (i - 0.50) / n,
    stop("Método no reconocido: ", metodo)
  )
}

calcular_aicc <- function(loglik, k, n) {
  aic <- -2 * loglik + 2 * k
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

ajustar_gumbel <- function(x) {
  s_ini  <- sqrt(6) * sd(x) / pi
  mu_ini <- mean(x) - EULER_MASCHERONI * s_ini
  nll <- function(par) {
    if (par[2] <= 0) return(Inf)
    -sum(evd::dgumbel(x, loc = par[1], scale = par[2], log = TRUE))
  }
  opt <- optim(c(mu_ini, s_ini), nll, method = "L-BFGS-B", lower = c(-Inf, 1e-6))
  list(mu = opt$par[1], sigma = opt$par[2], loglik = -opt$value)
}

ajustar_lognormal <- function(x) {
  fit <- fitdist(x, "lnorm")
  list(mu = fit$estimate["meanlog"], sigma = fit$estimate["sdlog"], loglik = fit$loglik)
}

ajustar_pearson3 <- function(x) {
  m1 <- mean(x); m2 <- var(x); cs <- moments::skewness(x)
  if (abs(cs) < 0.01) cs <- sign(cs) * 0.01
  shape0 <- (2/cs)^2; scale0 <- sqrt(m2/shape0)
  loc0   <- min(m1 - shape0*scale0, min(x) - 1e-4)
  nll <- function(par) {
    sh <- exp(par[1]); sc <- exp(par[2]); lo <- par[3]
    if (any(x <= lo)) return(1e12)
    -sum(dgamma(x - lo, shape = sh, scale = sc, log = TRUE))
  }
  opt <- optim(c(log(shape0), log(scale0), loc0), nll,
               method = "L-BFGS-B",
               lower  = c(-10, -10, -Inf),
               upper  = c(10,  10,  min(x) - 1e-8),
               control = list(maxit = 2000, factr = 1e7))
  if (opt$convergence != 0)
    warning("PE3: convergencia no garantizada (código ", opt$convergence, ")")
  par <- c(exp(opt$par[1]), exp(opt$par[2]), opt$par[3])
  list(shape = par[1], scale = par[2], loc = par[3],
       loglik = -opt$value, convergido = opt$convergence == 0)
}

ajustar_gamlss_seguro <- function(formula_mu, formula_sigma = ~1,
                                  familia, datos, nombre) {
  msgs <- character(0)
  m <- withCallingHandlers(
    tryCatch(
      gamlss(formula_mu, sigma.formula = formula_sigma,
             family = familia, data = datos, trace = FALSE),
      error = function(e) {
        cat(sprintf("  \u2717 %s: error \u2014 %s\n", nombre, e$message))
        NULL
      }
    ),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  if (!is.null(m) && length(msgs) > 0) {
    attr(m, "gamlss_warnings") <- msgs
    cat(sprintf("  \u26A0 %s: %d advertencia(s) de convergencia\n", nombre, length(msgs)))
  }
  m
}

bootstrap_gev_quantiles <- function(params, n_data, Tr_vec, B = 1000, alpha = 0.10) {
  set.seed(SEMILLA_GLOBAL)
  p_exc_vec <- 1 / Tr_vec
  boot_q <- matrix(NA_real_, B, length(Tr_vec))
  for (b in 1:B) {
    x_b <- revd(n_data, loc = params["location"],
                scale = params["scale"], shape = params["shape"])
    fit_b <- tryCatch(
      fevd(x_b, type = "GEV", method = "MLE", verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit_b))
      boot_q[b, ] <- qevd(1 - p_exc_vec,
                           loc   = fit_b$results$par["location"],
                           scale = fit_b$results$par["scale"],
                           shape = fit_b$results$par["shape"])
  }
  list(lower = apply(boot_q, 2, quantile, alpha/2,     na.rm = TRUE),
       upper = apply(boot_q, 2, quantile, 1-alpha/2,   na.rm = TRUE),
       B = B, alpha = alpha)
}

bootstrap_gamlss_quantiles <- function(modelo, datos, anios_pred, Tr_vec, B = 500) {
  set.seed(SEMILLA_GLOBAL)
  fam            <- modelo$family[1]
  n_obs          <- nrow(datos)
  p_exc_vec      <- 1 / Tr_vec
  media_anio_orig <- mean(datos$Anio)
  f_mu           <- modelo$mu.formula
  f_sigma        <- modelo$sigma.formula
  boot_mat <- array(NA_real_, dim = c(B, length(anios_pred), length(Tr_vec)))
  for (b in 1:B) {
    idx_b <- sample(n_obs, n_obs, replace = TRUE)
    db    <- datos[idx_b, ]
    m_b   <- tryCatch(
      suppressWarnings(gamlss(f_mu, sigma.formula = f_sigma,
                              family = fam, data = db, trace = FALSE)),
      error = function(e) NULL
    )
    if (is.null(m_b)) next
    mu_b_fit    <- tryCatch(as.numeric(fitted(m_b, what = "mu")),    error = function(e) NULL)
    sigma_b_fit <- tryCatch(as.numeric(fitted(m_b, what = "sigma")), error = function(e) NULL)
    if (is.null(mu_b_fit)) next
    for (i in seq_along(anios_pred)) {
      idx_m <- which(abs(db$Anio_centrado - (anios_pred[i] - media_anio_orig)) < 1e-8)
      if (length(idx_m) > 0) {
        mu_b <- mu_b_fit[idx_m[1]]; sigma_b <- sigma_b_fit[idx_m[1]]
      } else {
        nd <- data.frame(Anio_centrado = anios_pred[i] - media_anio_orig)
        mu_b <- tryCatch(
          as.numeric(predict(m_b, what = "mu",    newdata = nd, type = "response", data = db))[1],
          error = function(e) NA_real_)
        sigma_b <- tryCatch(
          as.numeric(predict(m_b, what = "sigma", newdata = nd, type = "response", data = db))[1],
          error = function(e) NA_real_)
      }
      if (anyNA(c(mu_b, sigma_b))) next
      for (j in seq_along(Tr_vec)) {
        boot_mat[b, i, j] <- switch(fam,
          "LOGNO" = qLOGNO(1-p_exc_vec[j], mu=mu_b, sigma=sigma_b),
          "GU"    = qGU   (1-p_exc_vec[j], mu=mu_b, sigma=sigma_b),
          "GA"    = qGA   (1-p_exc_vec[j], mu=mu_b, sigma=sigma_b),
          "NO"    = qNO   (1-p_exc_vec[j], mu=mu_b, sigma=sigma_b),
          NA_real_)
      }
    }
  }
  list(lower = apply(boot_mat, c(2,3), quantile, 0.025, na.rm = TRUE),
       upper = apply(boot_mat, c(2,3), quantile, 0.975, na.rm = TRUE),
       B = B)
}

# 0.3 Constantes ---------------------------------------------------------------
SEMILLA_GLOBAL   <- 2024L
EULER_MASCHERONI <- 0.5772156649
CUNNANE_A        <- 0.40
CUNNANE_B        <- 0.20
GRINGORTEN_A     <- 0.44
GRINGORTEN_B     <- 0.12
METODO_TRAZADO   <- "cunnane"
W_VENTANA_MIN    <- 10L
N_SIM_WAVELET    <- 1000L
B_BOOT_EST       <- 1000L
B_BOOT_GAMLSS    <- 500L
ALPHA_BOOT       <- 0.10
DELTA_AIC_UMBRAL <- 2.0
set.seed(SEMILLA_GLOBAL)

nivel_significancia <- 0.05
meses_espanol <- c("Enero","Febrero","Marzo","Abril","Mayo","Junio",
                   "Julio","Agosto","Septiembre","Octubre","Noviembre","Diciembre")
Tr_diseno <- c(2, 2.33, 5, 10, 15, 20, 25, 50, 100, 200, 500)
p_exc     <- 1 / Tr_diseno

CFG_GRAF <- list(dpi = 300L, ancho_s = 10, ancho_d = 14, alto_p = 12)

# 0.4 Rutas --------------------------------------------------------------------
ruta_archivo      <- "D:/R/Version_Final_CFE_Colombia/Qmax.xlsx"
directorio_salida <- "D:/R/Version_Final_CFE_Colombia/outputs"
dir.create(directorio_salida, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. CARGA Y ESTRUCTURACIÓN DE DATOS
# ==============================================================================
cat("\n", strrep("=", 70), "\n")
cat("   ANÁLISIS DE FRECUENCIAS DE CAUDALES MÁXIMOS — CFE Colombia\n")
cat("   Río Guachicono, Cauca — Colombia\n")
cat("   Autor : Mauricio Javier Victoria Niño\n")
cat("   ORCID : 0009-0003-4328-5691\n")
cat(strrep("=", 70), "\n\n")

datos_raw <- read_excel(ruta_archivo, sheet = 1)
if (length(names(datos_raw)) == 13)
  names(datos_raw) <- c("Anio", meses_espanol)

datos_mensuales <- datos_raw |>
  pivot_longer(cols = -Anio, names_to = "Mes", values_to = "Qmax") |>
  mutate(
    Anio    = as.integer(Anio),
    Mes_num = match(Mes, meses_espanol),
    Fecha   = as.Date(paste(Anio, Mes_num, "01", sep = "-")),
    Qmax    = as.numeric(Qmax)
  ) |>
  filter(!is.na(Qmax)) |>
  arrange(Fecha)

datos_anuales <- datos_mensuales |>
  group_by(Anio) |>
  summarise(
    Qmax_anual   = max(Qmax),
    Qmedia_anual = mean(Qmax),
    Qsd_anual    = sd(Qmax),
    Qcv_anual    = sd(Qmax) / mean(Qmax),
    .groups = "drop"
  ) |>
  arrange(Anio)

anios             <- datos_anuales$Anio
n_anios           <- length(anios)
serie_media_anual <- datos_anuales$Qmedia_anual
serie_max_anual   <- datos_anuales$Qmax_anual

cat("Datos cargados:\n")
cat("  Período    :", min(anios), "-", max(anios), "(", n_anios, "años)\n")
cat("  Mensuales  :", nrow(datos_mensuales), "registros\n")
cat("  Qmax anual :", round(min(serie_max_anual),2), "-",
    round(max(serie_max_anual),2), "m³/s | Media:",
    round(mean(serie_max_anual),2), "m³/s\n\n")


# ==============================================================================
# 2. MÓDULO 1: DETECCIÓN DE NO ESTACIONARIEDAD
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("   MÓDULO 1: DETECCIÓN DE NO ESTACIONARIEDAD\n")
cat(strrep("=", 70), "\n\n")

# 2.1 Análisis exploratorio visual --------------------------------------------
cat("  2.1 Panel exploratorio...\n")

datos_mensuales <- datos_mensuales |>
  mutate(
    Mes_factor = factor(Mes_num, 1:12, substr(meses_espanol, 1, 3)),
    Decada     = paste0(floor(Anio/10)*10, "s"),
    Epoca      = case_when(
      Mes_num %in% c(12,1,2)  ~ "Dic-Feb (Seca mayor)",
      Mes_num %in% c(3,4,5)   ~ "Mar-May (Transición)",
      Mes_num %in% c(6,7,8)   ~ "Jun-Ago (Seca menor)",
      Mes_num %in% c(9,10,11) ~ "Sep-Nov (Húmeda)"
    )
  )

loess_media <- loess(Qmedia_anual ~ Anio, data = datos_anuales, span = 0.4)
loess_sd    <- loess(Qsd_anual    ~ Anio, data = datos_anuales, span = 0.4)
datos_anuales$Loess_media <- predict(loess_media)
datos_anuales$Loess_sd    <- predict(loess_sd)

g1 <- ggplot(datos_mensuales, aes(Fecha, Qmax)) +
  geom_line(color = "#bdbdbd", linewidth = 0.2, alpha = 0.6) +
  geom_line(data = datos_anuales,
            aes(as.Date(paste(Anio,7,"01",sep="-")), Loess_media),
            color = "#2166ac", linewidth = 1.2) +
  geom_ribbon(data = datos_anuales,
              aes(as.Date(paste(Anio,7,"01",sep="-")),
                  ymin = Loess_media - Loess_sd,
                  ymax = Loess_media + Loess_sd),
              fill = "#2166ac", alpha = 0.1, inherit.aes = FALSE) +
  labs(title = "Serie temporal completa",
       subtitle = "Azul = tendencia LOESS | Banda = ±1 DE",
       y = "Caudal máximo mensual (m³/s)", x = "") +
  theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold"))

g2 <- ggplot(datos_mensuales, aes(Mes_factor, Qmax, fill = Epoca)) +
  geom_boxplot(outlier.size = 0.8, width = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
  scale_fill_manual(values = c("#fc8d59","#ffffbf","#91bfdb","#2b83ba")) +
  labs(title = "Patrón estacional", y = "Caudal máximo mensual (m³/s)", x = "") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

g3 <- ggplot(datos_mensuales, aes(Decada, Qmax, fill = Decada)) +
  geom_boxplot(outlier.size = 0.8, width = 0.6, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
  scale_fill_brewer(palette = "RdYlBu") +
  labs(title = "Evolución decadal", y = "Caudal máximo mensual (m³/s)", x = "") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

g4 <- ggplot(datos_anuales, aes(Anio)) +
  geom_ribbon(aes(ymin = Loess_media - Loess_sd, ymax = Loess_media + Loess_sd),
              fill = "#2166ac", alpha = 0.1) +
  geom_point(aes(y = Qmedia_anual), color = "#2c7bb6", size = 2) +
  geom_line(aes(y = Loess_media), color = "#2166ac", linewidth = 1.2) +
  geom_smooth(aes(y = Qmedia_anual), method = "lm", se = TRUE,
              color = "#d7191c", linewidth = 0.7, linetype = "dashed", alpha = 0.2) +
  labs(title = "Tendencia en caudal medio anual",
       subtitle = "Azul = LOESS | Rojo = Regresión lineal",
       y = "Caudal medio anual (m³/s)") +
  theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold"))

panel_exp <- grid.arrange(g1, g2, g3, g4, ncol = 2,
  top = textGrob("ANÁLISIS EXPLORATORIO — RÍO GUACHICONO",
                 gp = gpar(fontface = "bold", fontsize = 15)))
ggsave(file.path(directorio_salida, "Figura1_Panel_Exploratorio.png"),
       panel_exp, width = CFG_GRAF$ancho_d, height = CFG_GRAF$alto_p, dpi = CFG_GRAF$dpi)
cat("    \u2713 Figura 1 guardada\n")

# 2.2 Tests de estacionariedad ------------------------------------------------
cat("  2.2 Tests de estacionariedad...\n")

lag_lb      <- min(10L, floor(n_anios/5))
lb_test     <- Box.test(serie_media_anual, lag = lag_lb, type = "Ljung-Box")
autocorr_sig <- lb_test$p.value < nivel_significancia
cat(sprintf("    Ljung-Box (lag=%d): p = %.4f %s\n", lag_lb, lb_test$p.value,
    ifelse(autocorr_sig, "— \u26A0 autocorrelación significativa",
                         "— independencia no rechazada")))

resultados_tests <- list()

mk_test  <- MannKendall(serie_media_anual)
sen_test <- sens.slope(serie_media_anual)
resultados_tests$MK <- list(
  nombre        = "Mann-Kendall",
  estadistico   = paste0("\u03c4 = ", round(mk_test$tau[1],3),
                          ", Sen = ", round(sen_test$estimates,3), " m\u00b3/s/a\u00f1o"),
  p_valor       = mk_test$sl[1],
  significativo = mk_test$sl[1] < nivel_significancia,
  interpretacion = ifelse(mk_test$sl[1] < nivel_significancia,
    paste("Tendencia", ifelse(mk_test$tau[1] > 0, "creciente", "decreciente")),
    "Sin tendencia monótona")
)

pt_test       <- pettitt.test(serie_media_anual)
anio_cambio   <- anios[pt_test$estimate]
media_antes   <- mean(serie_media_anual[anios <  anio_cambio])
media_despues <- mean(serie_media_anual[anios >= anio_cambio])
cambio_pct    <- round((media_despues - media_antes) / media_antes * 100, 1)
resultados_tests$Pettitt <- list(
  nombre        = "Mann-Whitney-Pettitt",
  estadistico   = paste0("Cambio en ", anio_cambio, " (\u0394 = ", cambio_pct, "%)"),
  p_valor       = pt_test$p.value,
  significativo = pt_test$p.value < nivel_significancia,
  interpretacion = ifelse(pt_test$p.value < nivel_significancia,
    paste("Cambio estructural abrupto en", anio_cambio), "Sin cambio estructural")
)

u_forward  <- calcular_u(serie_media_anual)
u_backward <- -rev(calcular_u(rev(serie_media_anual)))
cruces     <- which(abs(u_forward - u_backward) < 0.8 &
                    sign(u_forward) != sign(u_backward))
sneyers_sig  <- any(abs(u_forward) > 1.96, na.rm=TRUE) ||
                any(abs(u_backward) > 1.96, na.rm=TRUE)
anio_sneyers <- if (length(cruces) > 0) {
  anios[cruces[which.min(abs(u_forward[cruces] - u_backward[cruces]))]]
} else { NA_integer_ }
resultados_tests$Sneyers <- list(
  nombre        = "Mann-Kendall Sneyers",
  estadistico   = ifelse(!is.na(anio_sneyers),
                         paste0("Cruce en ", anio_sneyers), "Sin cruce detectado"),
  p_valor       = NA,
  significativo = sneyers_sig,
  interpretacion = ifelse(sneyers_sig, "Cambio progresivo detectado", "Serie homogénea")
)

w_mw            <- max(W_VENTANA_MIN, round(n_anios * 0.25))
pvals_mw        <- sapply(1:(n_anios-w_mw+1),
                          function(i) MannKendall(serie_media_anual[i:(i+w_mw-1)])$sl[1])
taus_mw         <- sapply(1:(n_anios-w_mw+1),
                          function(i) MannKendall(serie_media_anual[i:(i+w_mw-1)])$tau[1])
prop_sig_mw     <- mean(pvals_mw < nivel_significancia)
cambios_signo_mw <- sum(diff(sign(taus_mw)) != 0)
estabilidad_mw  <- ifelse(prop_sig_mw > 0.6 & cambios_signo_mw == 0, "Alta",
                   ifelse(prop_sig_mw > 0.4, "Media", "Baja"))
resultados_tests$MWMK <- list(
  nombre        = paste0("Moving-Window MK (w=", w_mw, ")"),
  estadistico   = paste0("Prop. significativas = ", round(prop_sig_mw*100,1), "%"),
  p_valor       = NA,
  significativo = prop_sig_mw > 0.4,
  interpretacion = paste("Estabilidad de tendencia:", estabilidad_mw)
)

modelo_lm  <- lm(serie_media_anual ~ anios)
white_test <- bptest(modelo_lm, ~ anios + I(anios^2), data = data.frame(anios = anios))
resultados_tests$White <- list(
  nombre        = "White (heterocedasticidad)",
  estadistico   = paste0("BP = ", round(white_test$statistic, 3)),
  p_valor       = white_test$p.value,
  significativo = white_test$p.value < nivel_significancia,
  interpretacion = ifelse(white_test$p.value < nivel_significancia,
    "Varianza no constante", "Varianza constante")
)

adf_orig <- adf.test(serie_media_anual, alternative = "stationary")
adf_res  <- adf.test(residuals(modelo_lm),  alternative = "stationary")
resultados_tests$ADF <- list(
  nombre        = "Dickey-Fuller Aumentado",
  estadistico   = paste0("ADF = ", round(adf_orig$statistic, 3)),
  p_valor       = adf_orig$p.value,
  significativo = adf_orig$p.value < nivel_significancia,
  interpretacion = ifelse(adf_orig$p.value < nivel_significancia,
    "Serie estacionaria (rechaza raíz unitaria)", "Posible raíz unitaria")
)

n_sig <- sum(sapply(resultados_tests, function(x) x$significativo))
cat("    Tests significativos:", n_sig, "/", length(resultados_tests), "\n")

g_sneyers <- ggplot(data.frame(Anio=anios, uf=u_forward, ub=u_backward), aes(Anio)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed",
             color = "#d7191c", linewidth = 0.7) +
  geom_line(aes(y = uf, color = "u(t) Forward"),   linewidth = 1) +
  geom_line(aes(y = ub, color = "u'(t) Backward"), linewidth = 1) +
  scale_color_manual(values = c("u(t) Forward"="#2166ac","u'(t) Backward"="#d7191c")) +
  labs(title    = "Test de Mann-Kendall Sneyers",
       subtitle = "Líneas rojas = límites \u00b11.96 (\u03b1=0.05)",
       y = "Estadístico u(t)", x = "A\u00f1o", color = "") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")
ggsave(file.path(directorio_salida, "Figura2_Test_Sneyers.png"),
       g_sneyers, width = CFG_GRAF$ancho_s, height = 5, dpi = CFG_GRAF$dpi)
cat("    \u2713 Figura 2 guardada\n")

# 2.3 Clasificación -----------------------------------------------------------
hay_tendencia          <- resultados_tests$MK$significativo
hay_cambio             <- resultados_tests$Pettitt$significativo
residuos_estacionarios <- adf_res$p.value < nivel_significancia

tests_formales <- c("MK","Pettitt","White","ADF")
n_sig_formal   <- sum(sapply(resultados_tests[tests_formales], function(x) x$significativo))
ruta_estacionaria <- n_sig_formal < 2

tipo_no_est <- if (n_sig <= 1) {
  "ESTACIONARIA"
} else if (hay_tendencia && residuos_estacionarios && !hay_cambio) {
  "TENDENCIA DETERMINISTA"
} else if (hay_cambio && !hay_tendencia) {
  "CAMBIO ESTRUCTURAL ABRUPTO"
} else if (hay_tendencia && hay_cambio) {
  "TENDENCIA + CAMBIO ESTRUCTURAL"
} else {
  "COMPLEJA / MIXTA"
}
confianza <- ifelse(tipo_no_est == "COMPLEJA / MIXTA", "MEDIA", "ALTA")

cat(sprintf("    Clasificación: %s (Confianza: %s)\n", tipo_no_est, confianza))
cat("    Ruta analítica:",
    ifelse(ruta_estacionaria, "A (Estacionaria)", "B (No estacionaria)"), "\n\n")
if (!ruta_estacionaria)
  warning("No estacionariedad detectada. Cuantiles estacionarios reportados solo como referencia.")

# 2.4 Tendencia por época hidrológica ----------------------------------------
cat("  2.3 Tendencia por época hidrológica...\n")

datos_epoca <- datos_mensuales |>
  group_by(Anio, Epoca) |>
  summarise(Qmedia = mean(Qmax), .groups = "drop")

resultados_epoca <- datos_epoca |>
  group_by(Epoca) |>
  summarise(
    tau        = round(MannKendall(Qmedia)$tau[1], 3),
    p_valor    = round(MannKendall(Qmedia)$sl[1], 4),
    pendiente  = round(sens.slope(Qmedia)$estimates, 4),
    cambio_pct = round((last(Qmedia)-first(Qmedia))/first(Qmedia)*100, 1),
    .groups = "drop"
  ) |>
  mutate(Signif = ifelse(p_valor < nivel_significancia, "\u2713", "\u2717"))

g_epoca <- ggplot(datos_epoca, aes(Anio, Qmedia, color = Epoca)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
  facet_wrap(~ Epoca, scales = "free_y", ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(title    = "Tendencia de caudal medio por época hidrológica",
       subtitle = paste0("Período: ", min(anios), "-", max(anios)),
       y = "Caudal medio (m\u00b3/s)", x = "A\u00f1o") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave(file.path(directorio_salida, "Figura3_Tendencia_Epoca.png"),
       g_epoca, width = 12, height = 8, dpi = CFG_GRAF$dpi)
cat("    \u2713 Figura 3 guardada\n")

# 2.5 Frecuencia de extremos en ventanas móviles ------------------------------
cat("  2.4 Frecuencia de extremos...\n")

q90 <- quantile(datos_mensuales$Qmax, 0.90)
q95 <- quantile(datos_mensuales$Qmax, 0.95)
q99 <- quantile(datos_mensuales$Qmax, 0.99)
w_ext <- W_VENTANA_MIN

extremos <- data.frame()
for (i in 1:(n_anios - w_ext + 1)) {
  sub <- datos_mensuales |> filter(Anio >= anios[i], Anio <= anios[i+w_ext-1])
  nm  <- nrow(sub)
  extremos <- rbind(extremos, data.frame(
    Anio_central = anios[i + floor(w_ext/2)],
    P90 = sum(sub$Qmax > q90)/nm*100,
    P95 = sum(sub$Qmax > q95)/nm*100,
    P99 = sum(sub$Qmax > q99)/nm*100
  ))
}

g_extremos <- ggplot(
    extremos |> pivot_longer(c(P90,P95,P99), names_to="Percentil", values_to="Frecuencia"),
    aes(Anio_central, Frecuencia, color = Percentil)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 0.5, linetype = "dashed") +
  scale_color_manual(values = c("P90"="#2c7bb6","P95"="#fdae61","P99"="#d7191c")) +
  labs(title    = paste0("Frecuencia de caudales extremos (ventana ", w_ext, " a\u00f1os)"),
       subtitle = paste0("Umbrales: P90=",round(q90,1)," P95=",round(q95,1),
                          " P99=",round(q99,1)," m\u00b3/s"),
       y = "Frecuencia de excedencia (%)", x = "A\u00f1o central de la ventana") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(directorio_salida, "Figura4_Extremos_Ventanas.png"),
       g_extremos, width = 12, height = 6, dpi = CFG_GRAF$dpi)
cat("    \u2713 Figura 4 guardada\n")

# 2.6 Análisis wavelet --------------------------------------------------------
cat("  2.5 Wavelet (n.sim =", N_SIM_WAVELET, ")...\n")

set.seed(SEMILLA_GLOBAL)
wt <- analyze.wavelet(
  data.frame(date=as.Date(paste(anios,"07","01",sep="-")), value=serie_media_anual),
  "value", loess.span = 0, dt = 1, dj = 1/20,
  lowerPeriod = 2, upperPeriod = floor(n_anios/3),
  make.pval = TRUE, n.sim = N_SIM_WAVELET, verbose = FALSE)

png(file.path(directorio_salida, "Figura5_Wavelet.png"),
    width = 14, height = 8, units = "in", res = CFG_GRAF$dpi)
wt.image(wt, "value", legend.params = list(lab = "Potencia wavelet"),
         periodlab = "Período (a\u00f1os)", timelab = "A\u00f1o",
         main = paste0("Análisis Wavelet — Caudal medio anual (",
                        min(anios),"-",max(anios),")"),
         plot.contour = TRUE, plot.ridge = TRUE, col.ridge = "white")
dev.off()
cat("    \u2713 Figura 5 guardada\n\n")


# ==============================================================================
# 3. MÓDULO 2: ANÁLISIS DE FRECUENCIAS
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("   MÓDULO 2: ANÁLISIS DE FRECUENCIAS DE CAUDALES MÁXIMOS\n")
cat(strrep("=", 70), "\n\n")
cat("  Ruta:", ifelse(ruta_estacionaria,"A (Estacionario)","B (No Estacionario)"), "\n\n")

# 3.1 Análisis estacionario ---------------------------------------------------
cat("  3.1 Ajuste de distribuciones estacionarias...\n")

ajustes <- list()

tryCatch({
  fit_gum <- ajustar_gumbel(serie_max_anual)
  ajustes$Gumbel <- list(
    nombre = "Gumbel (GUM)",
    parametros = paste0("\u03bc=",round(fit_gum$mu,3),", \u03c3=",round(fit_gum$sigma,3)),
    loglik = fit_gum$loglik, n_param = 2,
    cuantiles = evd::qgumbel(1-p_exc, loc=fit_gum$mu, scale=fit_gum$sigma))
}, error = function(e) cat("    \u26A0 Error Gumbel:", e$message, "\n"))

tryCatch({
  fit_lno <- ajustar_lognormal(serie_max_anual)
  ajustes$LogNormal <- list(
    nombre = "Log-Normal (LNO)",
    parametros = paste0("\u03bclog=",round(fit_lno$mu,3),", \u03c3log=",round(fit_lno$sigma,3)),
    loglik = fit_lno$loglik, n_param = 2,
    cuantiles = qlnorm(1-p_exc, meanlog=fit_lno$mu, sdlog=fit_lno$sigma))
}, error = function(e) cat("    \u26A0 Error Log-Normal:", e$message, "\n"))

tryCatch({
  fit_pe3 <- ajustar_pearson3(serie_max_anual)
  ajustes$PearsonIII <- list(
    nombre = "Pearson III (PE3)",
    parametros = paste0("\u03b1=",round(fit_pe3$shape,3),
                         ", \u03b2=",round(fit_pe3$scale,3),
                         ", \u03be=",round(fit_pe3$loc,3)),
    loglik = fit_pe3$loglik, n_param = 3,
    cuantiles = qgamma(1-p_exc, shape=fit_pe3$shape, scale=fit_pe3$scale) + fit_pe3$loc)
}, error = function(e) cat("    \u26A0 Error Pearson III:", e$message, "\n"))

tryCatch({
  mu_ini_lp3    <- mean(log(serie_max_anual))
  sigma_ini_lp3 <- sd(log(serie_max_anual))
  fit_lp3 <- gamlss(serie_max_anual ~ 1, family = LOGNO,
                    mu.start = mu_ini_lp3, sigma.start = sigma_ini_lp3, trace = FALSE)
  ajustes$LogPearsonIII <- list(
    nombre = "Log-Pearson III (LP3)",
    parametros = paste0("\u03bc=",round(fit_lp3$mu.coefficients,3),
                         ", \u03c3=",round(exp(fit_lp3$sigma.coefficients),3)),
    loglik = as.numeric(logLik(fit_lp3)), n_param = 2,
    cuantiles = qLOGNO(1-p_exc, mu=fit_lp3$mu.coefficients,
                        sigma=fit_lp3$sigma.coefficients))
}, error = function(e) cat("    \u26A0 Error LP3:", e$message, "\n"))

tryCatch({
  fit_gev <- fevd(serie_max_anual, type = "GEV", method = "MLE", verbose = FALSE)
  pgev <- fit_gev$results$par
  ajustes$GEV <- list(
    nombre = "GEV",
    parametros = paste0("\u03be=",round(pgev["location"],3),
                         ", \u03c3=",round(pgev["scale"],3),
                         ", \u03ba=",round(pgev["shape"],3)),
    loglik = -fit_gev$results$value, n_param = 3,
    cuantiles = qevd(1-p_exc, loc=pgev["location"], scale=pgev["scale"], shape=pgev["shape"]))
}, error = function(e) cat("    \u26A0 Error GEV:", e$message, "\n"))

tryCatch({
  fit_nor <- fitdist(serie_max_anual, "norm")
  ajustes$Normal <- list(
    nombre = "Normal (NOR)",
    parametros = paste0("\u03bc=",round(fit_nor$estimate["mean"],3),
                         ", \u03c3=",round(fit_nor$estimate["sd"],3)),
    loglik = fit_nor$loglik, n_param = 2,
    cuantiles = qnorm(1-p_exc, mean=fit_nor$estimate["mean"], sd=fit_nor$estimate["sd"]))
}, error = function(e) cat("    \u26A0 Error Normal:", e$message, "\n"))

# AIC, AICc, BIC
for (nm in names(ajustes)) {
  a <- ajustes[[nm]]
  a$AIC  <- round(-2*a$loglik + 2*a$n_param, 2)
  a$AICc <- round(calcular_aicc(a$loglik, a$n_param, n_anios), 2)
  a$BIC  <- round(-2*a$loglik + a$n_param*log(n_anios), 2)
  ajustes[[nm]] <- a
}

# Filliben en n puntos (posición Cunnane)
filliben_est <- data.frame(Modelo=character(), r=numeric(), r2=numeric(),
                            stringsAsFactors=FALSE)
for (nm in names(ajustes)) {
  a   <- ajustes[[nm]]
  p_i <- posicion_trazado(n_anios, METODO_TRAZADO)
  q_teo <- tryCatch(switch(nm,
    Gumbel        = evd::qgumbel(p_i, loc=fit_gum$mu, scale=fit_gum$sigma),
    LogNormal     = qlnorm(p_i, meanlog=fit_lno$mu, sdlog=fit_lno$sigma),
    PearsonIII    = qgamma(p_i, shape=fit_pe3$shape, scale=fit_pe3$scale)+fit_pe3$loc,
    LogPearsonIII = qLOGNO(p_i, mu=fit_lp3$mu.coefficients, sigma=fit_lp3$sigma.coefficients),
    GEV           = qevd(p_i, loc=pgev["location"], scale=pgev["scale"], shape=pgev["shape"]),
    Normal        = qnorm(p_i, mean=fit_nor$estimate["mean"], sd=fit_nor$estimate["sd"]),
    rep(NA_real_, n_anios)), error = function(e) rep(NA_real_, n_anios))
  r_f <- suppressWarnings(cor(sort(serie_max_anual), q_teo, use="complete.obs"))
  filliben_est <- rbind(filliben_est,
    data.frame(Modelo=a$nombre, r=round(r_f,4), r2=round(r_f^2,4)))
}

aicc_vals  <- sapply(ajustes, function(x) x$AICc)
mejor_aicc <- names(ajustes)[which.min(aicc_vals)]
delta_aicc <- round(aicc_vals - min(aicc_vals), 2)
cat("    Mejor modelo (AICc):", ajustes[[mejor_aicc]]$nombre, "\n")
cat("    Modelos competitivos (\u0394AICc < 2):",
    paste(names(delta_aicc[delta_aicc < DELTA_AIC_UMBRAL]), collapse=", "), "\n")

# IC bootstrap GEV
ic_gev <- NULL
if ("GEV" %in% names(ajustes)) {
  cat("    Calculando IC bootstrap GEV (B =", B_BOOT_EST, ")...\n")
  ic_gev <- bootstrap_gev_quantiles(pgev, n_anios, Tr_diseno, B=B_BOOT_EST, alpha=ALPHA_BOOT)
  cat("    \u2713 IC GEV listo\n")
}

# Gráfico ajuste estacionario
x_seq     <- seq(min(serie_max_anual), max(serie_max_anual), length.out = 300)
datos_emp <- data.frame(x=sort(serie_max_anual), y=posicion_trazado(n_anios, METODO_TRAZADO))
curvas    <- data.frame(x = x_seq)
if ("Gumbel"    %in% names(ajustes))
  curvas$Gumbel    <- evd::pgumbel(x_seq, loc=fit_gum$mu, scale=fit_gum$sigma)
if ("LogNormal" %in% names(ajustes))
  curvas$LogNormal <- plnorm(x_seq, meanlog=fit_lno$mu, sdlog=fit_lno$sigma)
if ("Normal"    %in% names(ajustes))
  curvas$Normal    <- pnorm(x_seq, mean=fit_nor$estimate["mean"], sd=fit_nor$estimate["sd"])
curvas_long <- curvas |>
  pivot_longer(-x, names_to="Distribucion", values_to="Prob") |> filter(!is.na(Prob))

g_ajuste <- ggplot() +
  geom_point(data=datos_emp, aes(x,y), color="#2c7bb6", size=2, alpha=0.7) +
  geom_line(data=curvas_long, aes(x, Prob, color=Distribucion), linewidth=1.2) +
  scale_color_manual(values=c("Gumbel"="#d7191c","LogNormal"="#fdae61","Normal"="#5e3c99")) +
  labs(title    = "Ajuste de Distribuciones Estacionarias — Máximos Anuales",
       subtitle = paste0("Posición de trazado: Cunnane | Mejor modelo (AICc): ",
                          ajustes[[mejor_aicc]]$nombre),
       x = "Caudal máximo anual (m\u00b3/s)", y = "Probabilidad acumulada") +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(directorio_salida, "Figura6_Ajuste_Estacionario.png"),
       g_ajuste, width=CFG_GRAF$ancho_s, height=7, dpi=CFG_GRAF$dpi)
cat("    \u2713 Figura 6 guardada\n\n")

# 3.2 Análisis no estacionario (GAMLSS) ---------------------------------------
modelos_gamlss <- NULL; cuantiles_ns <- NULL; ic_ns <- NULL
mejor_modelo   <- NULL; validacion_ns <- NULL; r_fill_ns <- NULL
ks_res         <- list(statistic=NA, p.value=NA)

if (!ruta_estacionaria) {

  cat("  3.2 Análisis no estacionario GAMLSS...\n")

  datos_gamlss <- data.frame(
    Qmax          = serie_max_anual,
    Anio          = anios,
    Anio_centrado = anios - mean(anios)
  )
  modelos_gamlss <- list()

  m0 <- ajustar_gamlss_seguro(Qmax ~ 1, ~1, NO, datos_gamlss, "M0_Normal_Estacionario")
  m1 <- ajustar_gamlss_seguro(Qmax ~ Anio_centrado, ~1, GU, datos_gamlss, "M1_Gumbel_mu_lineal")
  m2 <- ajustar_gamlss_seguro(Qmax ~ pb(Anio_centrado), ~1, GU, datos_gamlss, "M2_Gumbel_mu_spline")
  m3 <- ajustar_gamlss_seguro(Qmax ~ pb(Anio_centrado), ~ Anio_centrado, GU, datos_gamlss, "M3_Gumbel_mu_sigma")
  m4 <- ajustar_gamlss_seguro(Qmax ~ pb(Anio_centrado), ~ pb(Anio_centrado), LOGNO, datos_gamlss, "M4_LogNormal_mu_sigma")
  m5 <- ajustar_gamlss_seguro(Qmax ~ pb(Anio_centrado), ~1, GA, datos_gamlss, "M5_Gamma_mu")

  if (!is.null(m0)) modelos_gamlss$M0_Normal_Estacionario <- m0
  if (!is.null(m1)) modelos_gamlss$M1_Gumbel_mu_lineal    <- m1
  if (!is.null(m2)) modelos_gamlss$M2_Gumbel_mu_spline    <- m2
  if (!is.null(m3)) modelos_gamlss$M3_Gumbel_mu_sigma     <- m3
  if (!is.null(m4)) modelos_gamlss$M4_LogNormal_mu_sigma  <- m4
  if (!is.null(m5)) modelos_gamlss$M5_Gamma_mu            <- m5

  if (length(modelos_gamlss) > 0) {

    aics         <- sapply(modelos_gamlss, AIC)
    mejor_modelo <- names(which.min(aics))
    modelo_final <- modelos_gamlss[[mejor_modelo]]
    fam          <- modelo_final$family[1]

    cat("    Mejor modelo GAMLSS:", mejor_modelo,
        "(AIC:", round(min(aics),2), ")\n")

    tryCatch({
      sm_info <- getSmo(modelo_final, what = "mu")
      if (!is.null(sm_info))
        cat(sprintf("    edf spline \u03bc: %.2f | df.fit: %.2f\n",
                    sm_info$edf, modelo_final$df.fit))
    }, error = function(e) NULL)

    # Función auxiliar de predicción (fitted para puntos en muestra,
    # predict.gamlss para puntos fuera de muestra)
    predecir_params <- function(mod, anio_c, datos_ref) {
      idx_match <- which(abs(datos_ref$Anio_centrado - anio_c) < 1e-8)
      if (length(idx_match) > 0) {
        list(mu    = as.numeric(fitted(mod, what="mu")   [idx_match[1]]),
             sigma = as.numeric(fitted(mod, what="sigma")[idx_match[1]]))
      } else {
        nd <- data.frame(Anio_centrado = anio_c)
        list(
          mu = tryCatch(
            as.numeric(predict(mod, what="mu",    newdata=nd,
                               type="response", data=datos_ref))[1],
            error = function(e) NA_real_),
          sigma = tryCatch(
            as.numeric(predict(mod, what="sigma", newdata=nd,
                               type="response", data=datos_ref))[1],
            error = function(e) NA_real_)
        )
      }
    }

    # Validación Cox-Snell con valores ajustados
    n_obs     <- length(serie_max_anual)
    mu_fit    <- as.numeric(fitted(modelo_final, what = "mu"))
    sigma_fit <- as.numeric(fitted(modelo_final, what = "sigma"))
    residuos_cs <- numeric(n_obs)
    for (i in 1:n_obs) {
      residuos_cs[i] <- switch(fam,
        "GU"    = pGU   (serie_max_anual[i], mu=mu_fit[i], sigma=sigma_fit[i]),
        "LOGNO" = pLOGNO(serie_max_anual[i], mu=mu_fit[i], sigma=sigma_fit[i]),
        "GA"    = pGA   (serie_max_anual[i], mu=mu_fit[i], sigma=sigma_fit[i]),
        "NO"    = pNO   (serie_max_anual[i], mu=mu_fit[i], sigma=sigma_fit[i]),
        NA_real_)
    }
    residuos_exp <- -log(1 - residuos_cs)
    ks_res <- tryCatch(ks.test(residuos_exp, "pexp", rate=1),
                       error = function(e) list(statistic=NA, p.value=NA))
    p_fill_res   <- posicion_trazado(n_obs, "cunnane")
    teoricos_exp <- qexp(p_fill_res, rate = 1)
    r_fill_ns    <- suppressWarnings(cor(sort(residuos_exp), teoricos_exp, use="complete.obs"))
    validacion_ns <- if (!is.na(ks_res$p.value) && ks_res$p.value > 0.05) "Adecuado" else "Cuestionable"
    cal_ns <- ifelse(r_fill_ns > 0.99, "Excelente",
               ifelse(r_fill_ns > 0.97, "Bueno",
                ifelse(r_fill_ns > 0.95, "Aceptable", "Cuestionable")))
    cat("    Validación: KS p =", round(ks_res$p.value,4),
        "| Filliben r =", round(r_fill_ns,4), "|", cal_ns, "\n\n")

    # Cuantiles no estacionarios
    anios_pred <- unique(sort(c(seq(min(anios), max(anios), by=5), min(anios), max(anios))))
    cuantiles_ns <- matrix(NA_real_, nrow=length(anios_pred), ncol=length(Tr_diseno),
                           dimnames=list(as.character(anios_pred), paste0("Tr",Tr_diseno)))
    for (i in seq_along(anios_pred)) {
      tryCatch({
        params <- predecir_params(modelo_final, anios_pred[i]-mean(anios), datos_gamlss)
        for (j in seq_along(Tr_diseno)) {
          cuantiles_ns[i,j] <- switch(fam,
            "GU"    = qGU   (1-p_exc[j], mu=params$mu, sigma=params$sigma),
            "LOGNO" = qLOGNO(1-p_exc[j], mu=params$mu, sigma=params$sigma),
            "GA"    = qGA   (1-p_exc[j], mu=params$mu, sigma=params$sigma),
            "NO"    = qNO   (1-p_exc[j], mu=params$mu, sigma=params$sigma),
            NA_real_)
        }
      }, error = function(e) NULL)
    }

    # IC bootstrap GAMLSS
    cat("    Calculando IC bootstrap GAMLSS (B =", B_BOOT_GAMLSS, ")...\n")
    ic_ns <- bootstrap_gamlss_quantiles(modelo_final, datos_gamlss, anios_pred,
                                         Tr_diseno, B=B_BOOT_GAMLSS)
    cat("    \u2713 IC GAMLSS listo\n")

    # Figura 7: Q-Q Cox-Snell
    g_qq <- ggplot(data.frame(Teorico=teoricos_exp, Observado=sort(residuos_exp)),
                   aes(Teorico, Observado)) +
      geom_point(color="#2c7bb6", size=2.5, alpha=0.7) +
      geom_abline(slope=1, intercept=0, color="#d7191c", linewidth=0.8, linetype="dashed") +
      annotate("text", x=max(teoricos_exp)*0.3, y=max(residuos_exp)*0.95,
               label=paste0("r = ",round(r_fill_ns,4),"\nKS p = ",round(ks_res$p.value,4)),
               hjust=0, size=3.5) +
      labs(title    = "Validación del Modelo No Estacionario — Q-Q Plot Cox-Snell",
           subtitle = paste0("Modelo: ",mejor_modelo," | Línea roja = Exp(1)"),
           x = "Cuantiles teóricos Exp(1)", y = "Residuos de Cox-Snell ordenados") +
      theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold"))
    ggsave(file.path(directorio_salida,"Figura7_Validacion_NoEstacionario.png"),
           g_qq, width=9, height=7, dpi=CFG_GRAF$dpi)
    cat("    \u2713 Figura 7 guardada\n")

    # Figura 8: Evolución temporal de cuantiles con IC
    idx_list <- list(); etq_list <- list()
    for (tr_val in c(2.33, 10, 25, 50, 100)) {
      idx <- encontrar_idx(tr_val, Tr_diseno)
      if (!is.na(idx)) { idx_list <- c(idx_list, idx); etq_list <- c(etq_list, paste0("Tr = ",tr_val," años")) }
    }
    ref_dist <- NULL; ref_nombre <- ""
    if ("GEV" %in% names(ajustes)) { ref_dist <- ajustes$GEV$cuantiles; ref_nombre <- "GEV" }
    else if ("LogNormal" %in% names(ajustes)) { ref_dist <- ajustes$LogNormal$cuantiles; ref_nombre <- "Log-Normal" }

    if (length(idx_list) > 0) {
      n_tr <- length(idx_list); n_pred <- length(anios_pred)
      cuantiles_plot <- data.frame(
        Anio=rep(anios_pred,n_tr), Caudal=numeric(n_pred*n_tr),
        IC_lower=numeric(n_pred*n_tr), IC_upper=numeric(n_pred*n_tr),
        Tr_label=rep(unlist(etq_list), each=n_pred), stringsAsFactors=FALSE)
      for (m in seq_len(n_tr)) {
        j <- idx_list[[m]]; idx_rows <- ((m-1)*n_pred+1):(m*n_pred)
        cuantiles_plot$Caudal  [idx_rows] <- cuantiles_ns[,j]
        cuantiles_plot$IC_lower[idx_rows] <- if (!is.null(ic_ns)) ic_ns$lower[,j] else NA
        cuantiles_plot$IC_upper[idx_rows] <- if (!is.null(ic_ns)) ic_ns$upper[,j] else NA
      }

      g_ns <- ggplot() +
        geom_ribbon(data=cuantiles_plot, aes(Anio, ymin=IC_lower, ymax=IC_upper, fill=Tr_label), alpha=0.15) +
        geom_line(data=cuantiles_plot, aes(Anio, Caudal, color=Tr_label), linewidth=1.2) +
        geom_point(data=datos_anuales, aes(Anio, Qmax_anual), color="grey50", size=2, alpha=0.5) +
        facet_wrap(~ Tr_label, scales="free_y", ncol=2) +
        scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1") +
        labs(title    = "Cuantiles No Estacionarios — Evolución Temporal",
             subtitle = paste0("Modelo: ",mejor_modelo," | Banda = IC 95% bootstrap"),
             y = "Caudal máximo anual (m\u00b3/s)", x = "A\u00f1o") +
        theme_minimal(base_size=11) + theme(legend.position="none", plot.title=element_text(face="bold"))

      if (!is.null(ref_dist)) {
        ref_est <- data.frame(Anio_start=min(anios_pred), Anio_end=max(anios_pred),
                              Caudal=sapply(idx_list, function(j) ref_dist[j]),
                              Tr_label=unlist(etq_list), stringsAsFactors=FALSE)
        g_ns <- g_ns +
          geom_segment(data=ref_est, aes(x=Anio_start, xend=Anio_end, y=Caudal, yend=Caudal, color=Tr_label),
                       linetype="dashed", linewidth=0.7, alpha=0.6) +
          labs(subtitle=paste0("Modelo: ",mejor_modelo," | Punteado = ",ref_nombre,
                                " estacionario | Banda = IC 95% bootstrap"))
      }
      ggsave(file.path(directorio_salida,"Figura8_Cuantiles_NoEstacionarios.png"),
             g_ns, width=CFG_GRAF$ancho_d, height=CFG_GRAF$alto_p, dpi=CFG_GRAF$dpi)
      cat("    \u2713 Figura 8 guardada\n\n")
    }
  }
}


# ==============================================================================
# 4. TABLA RESUMEN FINAL (consola)
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("   RESUMEN — CAUDALES POR PERÍODO DE RETORNO (m³/s)\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("   %-6s | %-13s | %-13s | %-13s\n",
    "Tr", "GEV (Estac)", "LogNormal (Estac)",
    if (!is.null(cuantiles_ns)) paste0("NoEst (",max(anios),")") else ""))
cat("   ", strrep("-", 55), "\n")
for (i in seq_along(Tr_diseno)) {
  gev_v <- if ("GEV"       %in% names(ajustes)) sprintf("%.2f", ajustes$GEV$cuantiles[i])       else "N/A"
  lno_v <- if ("LogNormal" %in% names(ajustes)) sprintf("%.2f", ajustes$LogNormal$cuantiles[i]) else "N/A"
  cat(sprintf("   %-6s | %-13s | %-13s",
      ifelse(Tr_diseno[i]==2.33,"2.33",as.character(Tr_diseno[i])), gev_v, lno_v))
  if (!is.null(cuantiles_ns)) {
    idx_actual <- which(rownames(cuantiles_ns) == as.character(max(anios)))
    cat(sprintf(" | %-13.2f", cuantiles_ns[idx_actual, i]))
  }
  cat("\n")
}
cat("   ", strrep("-", 55), "\n\n")


# ==============================================================================
# 5. INFORME FINAL
# ==============================================================================
cat("Generando informe final...\n")

inf <- c()
inf <- c(inf, strrep("\u2550", 80))
inf <- c(inf, "   INFORME FINAL \u2014 ANÁLISIS DE FRECUENCIAS DE CAUDALES MÁXIMOS")
inf <- c(inf, "   CFE Colombia v1.0.0")
inf <- c(inf, strrep("\u2550", 80))
inf <- c(inf, "")
inf <- c(inf, paste("Fecha del análisis:", as.character(Sys.Date())))
inf <- c(inf, "Estación  : Río Guachicono, Cauca \u2014 Colombia")
inf <- c(inf, "Variable  : Caudal máximo mensual y anual (m\u00b3/s)")
inf <- c(inf, paste("Período   :", min(anios), "-", max(anios), "(", n_anios, "años)"))
inf <- c(inf, paste("Registros mensuales:", nrow(datos_mensuales)))
inf <- c(inf, "")
inf <- c(inf, "Autor     : Mauricio Javier Victoria Niño")
inf <- c(inf, "Filiación : Investigador independiente \u00b7 Cali, Colombia")
inf <- c(inf, "Contacto  : hidratecsa@gmail.com")
inf <- c(inf, "ORCID     : 0009-0003-4328-5691")
inf <- c(inf, "")

# Módulo 1
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, "MÓDULO 1: DETECCIÓN Y CARACTERIZACIÓN DE NO ESTACIONARIEDAD")
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, "")
inf <- c(inf, "1.1 ESTADÍSTICAS DE LA SERIE")
inf <- c(inf, "")
inf <- c(inf, "   Caudal máximo mensual (m\u00b3/s):")
inf <- c(inf, sprintf("   \u2022 Mínimo: %.2f | Máximo: %.2f | Media: %.2f | Mediana: %.2f",
  min(datos_mensuales$Qmax), max(datos_mensuales$Qmax),
  mean(datos_mensuales$Qmax), median(datos_mensuales$Qmax)))
inf <- c(inf, sprintf("   \u2022 Desviación estándar: %.2f | CV: %.3f | Asimetría: %.3f",
  sd(datos_mensuales$Qmax), sd(datos_mensuales$Qmax)/mean(datos_mensuales$Qmax),
  skewness(datos_mensuales$Qmax)))
inf <- c(inf, "")
inf <- c(inf, "   Caudal máximo anual (m\u00b3/s):")
inf <- c(inf, sprintf("   \u2022 Mínimo: %.2f (%d) | Máximo: %.2f (%d) | Media: %.2f",
  min(serie_max_anual), anios[which.min(serie_max_anual)],
  max(serie_max_anual), anios[which.max(serie_max_anual)], mean(serie_max_anual)))
inf <- c(inf, sprintf("   \u2022 CV: %.3f | Asimetría: %.3f",
  sd(serie_max_anual)/mean(serie_max_anual), skewness(serie_max_anual)))
inf <- c(inf, "")

inf <- c(inf, "1.2 TEST PREVIO DE AUTOCORRELACIÓN")
inf <- c(inf, "")
inf <- c(inf, sprintf("   Ljung-Box (lag=%d): estadístico = %.4f, p-valor = %.4f  %s",
  lag_lb, lb_test$statistic, lb_test$p.value,
  ifelse(autocorr_sig, "\u26A0 Autocorrelación significativa", "\u2713 Independencia no rechazada")))
inf <- c(inf, "")

inf <- c(inf, "1.3 BATERÍA DE TESTS DE ESTACIONARIEDAD (\u03b1 = 0.05)")
inf <- c(inf, "")
inf <- c(inf, "   Serie analizada: Caudal medio anual")
inf <- c(inf, "")
for (r in resultados_tests) {
  inf <- c(inf, paste("   \u250C\u2500", r$nombre))
  inf <- c(inf, paste("   \u2502 Estadístico:", r$estadistico))
  if (!is.na(r$p_valor)) {
    p_str <- sprintf("%.4f", r$p_valor)
    if (r$p_valor < 0.01) p_str <- paste0(p_str," (evidencia MUY FUERTE)")
    else if (r$p_valor < 0.05) p_str <- paste0(p_str," (evidencia SIGNIFICATIVA)")
    inf <- c(inf, paste("   \u2502 p-valor:", p_str))
  }
  inf <- c(inf, paste("   \u2502 Resultado:", ifelse(r$significativo, "\u2713 SIGNIFICATIVO", "\u2717 NO SIGNIFICATIVO")))
  inf <- c(inf, paste("   \u2502 Interpretación:", r$interpretacion))
  inf <- c(inf, "   \u2514\u2500"); inf <- c(inf, "")
}
inf <- c(inf, sprintf("   RESUMEN: %d de %d tests significativos", n_sig, length(resultados_tests)))
inf <- c(inf, sprintf("   Tests formales (MK/Pettitt/White/ADF): %d / 4 significativos", n_sig_formal))
inf <- c(inf, "   Criterio de routing: \u2265 2 de 4 tests formales \u2192 No estacionaria")
inf <- c(inf, "   (Villarini et al. 2009, Adv. Water Resour. 32:1255-1266)")
inf <- c(inf, ""); inf <- c(inf, paste("   CLASIFICACIÓN FINAL:", tipo_no_est))
inf <- c(inf, paste("   Nivel de confianza:", confianza))
inf <- c(inf, "")

if (hay_tendencia) {
  inf <- c(inf, "   Detalle Mann-Kendall:")
  inf <- c(inf, sprintf("   \u2022 \u03c4 de Kendall: %.3f", mk_test$tau[1]))
  inf <- c(inf, sprintf("   \u2022 Pendiente de Sen: %.4f m\u00b3/s por año", sen_test$estimates))
  inf <- c(inf, sprintf("   \u2022 Cambio total estimado: %.2f m\u00b3/s en %d años", sen_test$estimates*n_anios, n_anios))
  inf <- c(inf, paste("   \u2022 Dirección:", ifelse(mk_test$tau[1]>0,"incremento","disminución")))
  inf <- c(inf, "")
}
if (hay_cambio) {
  inf <- c(inf, "   Detalle Pettitt:")
  inf <- c(inf, sprintf("   \u2022 Punto de cambio: %d", anio_cambio))
  inf <- c(inf, sprintf("   \u2022 Media antes: %.2f m\u00b3/s", media_antes))
  inf <- c(inf, sprintf("   \u2022 Media después: %.2f m\u00b3/s", media_despues))
  inf <- c(inf, sprintf("   \u2022 Magnitud del cambio: %.1f%%", cambio_pct))
  inf <- c(inf, "")
}

inf <- c(inf, "1.4 ANÁLISIS DE TENDENCIA POR ÉPOCA HIDROLÓGICA")
inf <- c(inf, "")
inf <- c(inf, sprintf("   %-25s | %8s | %8s | %10s | %8s", "Época","\u03c4","p-valor","Pendiente","Cambio%"))
inf <- c(inf, paste0("   ", strrep("-", 70)))
for (i in 1:nrow(resultados_epoca)) {
  inf <- c(inf, sprintf("   %-25s | %8.3f | %8.4f | %10.4f | %7.1f%%   %s",
    resultados_epoca$Epoca[i], resultados_epoca$tau[i], resultados_epoca$p_valor[i],
    resultados_epoca$pendiente[i], resultados_epoca$cambio_pct[i], resultados_epoca$Signif[i]))
}
inf <- c(inf, "")
inf <- c(inf, "1.5 FRECUENCIA DE EVENTOS EXTREMOS")
inf <- c(inf, "")
inf <- c(inf, paste("   Ventana móvil:", w_ext, "años"))
inf <- c(inf, sprintf("   Umbrales globales: P90=%.1f | P95=%.1f | P99=%.1f m\u00b3/s", q90,q95,q99))
inf <- c(inf, "")
inf <- c(inf, "1.6 ANÁLISIS WAVELET MULTI-ESCALA")
inf <- c(inf, paste0("   Wavelet madre Morlet | n.sim = ", N_SIM_WAVELET, " (Monte Carlo)"))
inf <- c(inf, "   (Ver Figura 5)")
inf <- c(inf, "")

# Módulo 2
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, "MÓDULO 2: ANÁLISIS DE FRECUENCIAS DE CAUDALES MÁXIMOS")
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, "")
inf <- c(inf, paste("   Ruta de análisis:", ifelse(ruta_estacionaria,"A (Estacionario)","B (No Estacionario)")))
inf <- c(inf, "")
inf <- c(inf, "2.1 ANÁLISIS ESTACIONARIO")
inf <- c(inf, "")
inf <- c(inf, "   Posición de trazado: Cunnane (propósito general)")
inf <- c(inf, paste("   Criterio de selección: AICc (corrección para n =", n_anios, ")"))
inf <- c(inf, "")
inf <- c(inf, sprintf("   %-22s | %8s | %8s | %8s | %8s | %8s",
                       "Distribución","AIC","AICc","BIC","Fill.r","Fill.r²"))
inf <- c(inf, paste0("   ", strrep("-", 75)))
for (nm in names(ajustes)) {
  a   <- ajustes[[nm]]
  fr  <- filliben_est$r [filliben_est$Modelo == a$nombre]
  fr2 <- filliben_est$r2[filliben_est$Modelo == a$nombre]
  inf <- c(inf, sprintf("   %-22s | %8.2f | %8.2f | %8.2f | %8.4f | %8.4f %s",
    a$nombre, a$AIC, a$AICc, a$BIC, fr, fr2, ifelse(nm==mejor_aicc,"<- mejor","")))
}
inf <- c(inf, "")
inf <- c(inf, paste("   Mejor modelo (AICc):", ajustes[[mejor_aicc]]$nombre))
inf <- c(inf, paste("   Parámetros:", ajustes[[mejor_aicc]]$parametros))
inf <- c(inf, sprintf("   Modelos con \u0394AICc < %.1f (competitivos): %s",
  DELTA_AIC_UMBRAL, paste(names(delta_aicc[delta_aicc < DELTA_AIC_UMBRAL]), collapse=", ")))
inf <- c(inf, "")
inf <- c(inf, "   Caudales estacionarios por período de retorno (m\u00b3/s):")
if (!is.null(ic_gev))
  inf <- c(inf, sprintf("   IC al %.0f%% bootstrap paramétrico GEV (B = %d réplicas)",
                         (1-ALPHA_BOOT)*100, ic_gev$B))
inf <- c(inf, sprintf("   %-6s | %-12s | %-12s | %-12s | %-14s | %-14s",
                       "Tr","LP3","GEV","Log-Normal","IC_inf GEV","IC_sup GEV"))
inf <- c(inf, paste0("   ", strrep("-", 80)))
for (i in seq_along(Tr_diseno)) {
  lp3_v <- if ("LogPearsonIII" %in% names(ajustes)) sprintf("%.2f", ajustes$LogPearsonIII$cuantiles[i]) else "N/A"
  gev_v <- if ("GEV"           %in% names(ajustes)) sprintf("%.2f", ajustes$GEV$cuantiles[i])           else "N/A"
  lno_v <- if ("LogNormal"     %in% names(ajustes)) sprintf("%.2f", ajustes$LogNormal$cuantiles[i])     else "N/A"
  ic_i  <- if (!is.null(ic_gev)) sprintf("%.2f", ic_gev$lower[i]) else "N/A"
  ic_s  <- if (!is.null(ic_gev)) sprintf("%.2f", ic_gev$upper[i]) else "N/A"
  inf <- c(inf, sprintf("   %-6s | %-12s | %-12s | %-12s | %-14s | %-14s",
    ifelse(Tr_diseno[i]==2.33,"2.33",as.character(Tr_diseno[i])), lp3_v, gev_v, lno_v, ic_i, ic_s))
}
inf <- c(inf, "")

if (!is.null(modelos_gamlss) && length(modelos_gamlss) > 0) {
  inf <- c(inf, "2.2 ANÁLISIS NO ESTACIONARIO (GAMLSS)")
  inf <- c(inf, "")
  inf <- c(inf, "   Modelos evaluados:")
  inf <- c(inf, sprintf("   %-30s | %8s | %8s | %s", "Modelo","AIC","BIC","Advertencias"))
  inf <- c(inf, paste0("   ", strrep("-", 65)))
  for (nm in names(modelos_gamlss)) {
    m    <- modelos_gamlss[[nm]]
    wstr <- if (!is.null(attr(m,"gamlss_warnings"))) paste0(length(attr(m,"gamlss_warnings"))," warn.") else "OK"
    inf <- c(inf, sprintf("   %-30s | %8.2f | %8.2f | %s", nm, AIC(m), BIC(m), wstr))
  }
  inf <- c(inf, "")
  inf <- c(inf, paste("   Mejor modelo:", mejor_modelo))
  inf <- c(inf, paste("   AIC:", round(min(aics),2)))
  inf <- c(inf, "")
  inf <- c(inf, "   Validación (Residuos de Cox-Snell):")
  inf <- c(inf, sprintf("   \u2022 Test KS: estadístico = %.4f, p-valor = %.4f", ks_res$statistic, ks_res$p.value))
  inf <- c(inf, sprintf("   \u2022 Filliben (n puntos, Cunnane): r = %.4f, r² = %.4f (%s)", r_fill_ns, r_fill_ns^2, cal_ns))
  inf <- c(inf, sprintf("   \u2022 Validación global: %s", validacion_ns))
  inf <- c(inf, "")
  inf <- c(inf, "   Cuantiles no estacionarios (m\u00b3/s) — estimación puntual:")
  cab <- sprintf("   %6s |", "Año")
  for (tr in Tr_diseno) cab <- paste0(cab, sprintf(" %7s |", paste0("Tr",tr)))
  inf <- c(inf, paste0("   ", strrep("-", 130)))
  inf <- c(inf, cab)
  inf <- c(inf, paste0("   ", strrep("-", 130)))
  for (i in 1:nrow(cuantiles_ns)) {
    lin <- sprintf("   %6d |", anios_pred[i])
    for (j in 1:ncol(cuantiles_ns)) lin <- paste0(lin, sprintf(" %7.2f |", cuantiles_ns[i,j]))
    inf <- c(inf, lin)
  }
  inf <- c(inf, paste0("   ", strrep("-", 130)))
  inf <- c(inf, "")

  idx_actual <- which(rownames(cuantiles_ns) == as.character(max(anios)))
  if (!is.null(ic_ns)) {
    inf <- c(inf, sprintf("   Intervalos de confianza 95%% bootstrap residuos GAMLSS (B = %d):", ic_ns$B))
    inf <- c(inf, "   (año más reciente de la serie)")
    inf <- c(inf, sprintf("   %-6s | %-12s | %-14s | %-14s", "Tr","Puntual","IC_inf 2.5%","IC_sup 97.5%"))
    inf <- c(inf, paste0("   ", strrep("-", 55)))
    for (i in seq_along(Tr_diseno)) {
      inf <- c(inf, sprintf("   %-6s | %-12.2f | %-14.2f | %-14.2f",
        ifelse(Tr_diseno[i]==2.33,"2.33",as.character(Tr_diseno[i])),
        cuantiles_ns[idx_actual,i], ic_ns$lower[idx_actual,i], ic_ns$upper[idx_actual,i]))
    }
    inf <- c(inf, "")
  }

  inf <- c(inf, "   Comparación con GEV estacionario (año más reciente):")
  inf <- c(inf, sprintf("   %-6s | %-15s | %-15s | %-15s", "Tr","No Estacionario","Estacionario","Diferencia"))
  inf <- c(inf, paste0("   ", strrep("-", 60)))
  for (i in seq_along(Tr_diseno)) {
    dif     <- cuantiles_ns[idx_actual,i] - ajustes$GEV$cuantiles[i]
    dif_pct <- dif / ajustes$GEV$cuantiles[i] * 100
    inf <- c(inf, sprintf("   %-6s | %-15.2f | %-15.2f | %+.2f (%+.1f%%)",
      ifelse(Tr_diseno[i]==2.33,"2.33",as.character(Tr_diseno[i])),
      cuantiles_ns[idx_actual,i], ajustes$GEV$cuantiles[i], dif, dif_pct))
  }
  inf <- c(inf, paste0("   ", strrep("-", 60))); inf <- c(inf, "")
}

# Notas metodológicas
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, "NOTAS METODOLÓGICAS")
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, "")
inf <- c(inf, "Módulo 1:")
inf <- c(inf, "\u2022 Test previo de autocorrelación: Ljung-Box")
inf <- c(inf, "\u2022 Batería de 6 tests: MK, Pettitt, Sneyers, MW-MK, White, ADF")
inf <- c(inf, "\u2022 \u03b1 = 0.05 para todos los tests de hipótesis")
inf <- c(inf, "\u2022 Criterio de no estacionariedad: \u2265 2/4 tests formales significativos")
inf <- c(inf, "\u2022 Análisis por época hidrológica: 4 períodos intra-anuales")
inf <- c(inf, paste0("\u2022 Wavelet madre Morlet | n.sim = ", N_SIM_WAVELET, " simulaciones Monte Carlo"))
inf <- c(inf, "")
inf <- c(inf, "Módulo 2:")
inf <- c(inf, "\u2022 Distribuciones: LP3, GEV, Gumbel, Log-Normal, Pearson III, Normal")
inf <- c(inf, "\u2022 Estimación por Máxima Verosimilitud (MLE)")
inf <- c(inf, "\u2022 Posición de trazado: Cunnane (uniforme en gráficos y Filliben)")
inf <- c(inf, "\u2022 Selección de modelo: AICc (Burnham & Anderson 2002)")
inf <- c(inf, paste0("\u2022 IC cuantiles GEV: bootstrap paramétrico, B = ",B_BOOT_EST,", IC ",round((1-ALPHA_BOOT)*100),"%"))
inf <- c(inf, "\u2022 GAMLSS: Modelos Aditivos Generalizados para Locación, Escala y Forma")
inf <- c(inf, paste0("\u2022 IC cuantiles GAMLSS: bootstrap de residuos, B = ",B_BOOT_GAMLSS,", IC 95%"))
inf <- c(inf, "\u2022 Validación no estacionaria: Cox-Snell (fitted) + Filliben + KS")
inf <- c(inf, "\u2022 Filliben calculado en n puntos de posición de trazado (no interpolado)")
inf <- c(inf, paste0("\u2022 Períodos de retorno: ", paste(Tr_diseno, collapse=", "), " años"))
inf <- c(inf, "")
inf <- c(inf, "Figuras:")
inf <- c(inf, paste0("\u2022 DPI: ", CFG_GRAF$dpi, " (JoH / HESS / WRR)"))
inf <- c(inf, "\u2022 Figura 1: Panel exploratorio 4-en-1")
inf <- c(inf, "\u2022 Figura 2: Test de Mann-Kendall Sneyers (secuencial)")
inf <- c(inf, "\u2022 Figura 3: Tendencia por época hidrológica")
inf <- c(inf, "\u2022 Figura 4: Frecuencia de extremos en ventanas móviles")
inf <- c(inf, "\u2022 Figura 5: Análisis wavelet multi-escala")
inf <- c(inf, "\u2022 Figura 6: Ajuste de distribuciones estacionarias (Cunnane)")
if (!is.null(cuantiles_ns)) {
  inf <- c(inf, "\u2022 Figura 7: Q-Q Cox-Snell (validación no estacionaria)")
  inf <- c(inf, "\u2022 Figura 8: Evolución temporal de cuantiles (IC 95% bootstrap)")
}

# Sesión R
inf <- c(inf, "")
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, "INFORMACIÓN DE SESIÓN R")
inf <- c(inf, strrep("\u2500", 80))
inf <- c(inf, capture.output(sessionInfo()))
inf <- c(inf, "")
inf <- c(inf, strrep("\u2550", 80))
inf <- c(inf, "   FIN DEL INFORME")
inf <- c(inf, strrep("\u2550", 80))

writeLines(inf, file.path(directorio_salida, "Informe_Final.txt"), useBytes = FALSE)
cat("\u2713 Informe_Final.txt guardado\n\n")


# ==============================================================================
# 6. VERIFICACIÓN DE ARCHIVOS GENERADOS
# ==============================================================================
cat(strrep("=", 70), "\n")
cat("   VERIFICACIÓN DE ARCHIVOS GENERADOS\n")
cat(strrep("=", 70), "\n\n")

archivos_esperados <- c(
  "Figura1_Panel_Exploratorio.png", "Figura2_Test_Sneyers.png",
  "Figura3_Tendencia_Epoca.png",    "Figura4_Extremos_Ventanas.png",
  "Figura5_Wavelet.png",            "Figura6_Ajuste_Estacionario.png",
  "Figura7_Validacion_NoEstacionario.png",
  "Figura8_Cuantiles_NoEstacionarios.png",
  "Informe_Final.txt"
)
for (archivo in archivos_esperados) {
  ruta_c <- file.path(directorio_salida, archivo)
  if (file.exists(ruta_c)) {
    cat(sprintf("  \u2713 %-48s (%.1f KB)\n", archivo, file.info(ruta_c)$size/1024))
  } else {
    cat(sprintf("  \u2717 %-48s NO ENCONTRADO\n", archivo))
  }
}

cat("\n", strrep("=", 70), "\n")
cat("   ANÁLISIS COMPLETADO EXITOSAMENTE\n")
cat("   Autor: Mauricio Javier Victoria Niño | ORCID: 0009-0003-4328-5691\n")
cat(strrep("=", 70), "\n")

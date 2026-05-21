# CFE-Colombia v1.0.0

## Análisis de Frecuencias de Caudales Máximos bajo No Estacionariedad

[![Licencia: CC BY 4.0](https://img.shields.io/badge/Licencia-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/deed.es)
[![R versión](https://img.shields.io/badge/R-%3E%3D4.3.1-blue.svg)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-10.XXXX%2Fengrxiv.XXXXXXXX-orange.svg)](https://doi.org/10.XXXX/engrxiv.XXXXXXXX)

**Autor:** Mauricio Javier Victoria Niño  
**Filiación:** Investigador independiente · Cali, Colombia  
**Contacto:** hidratecsa@gmail.com  
**ORCID:** [0009-0003-4328-5691](https://orcid.org/0009-0003-4328-5691)  
**Versión:** 1.0.0 · Publicado el 2026-05-20  

---

## Descripción general

**CFE-Colombia v1.0.0** es un marco computacional de código abierto en R para el análisis de frecuencias de caudales máximos bajo no estacionariedad. Fue desarrollado y validado en la cuenca alta del río Guachicono (estación IDEAM 52027010, Cauca, Colombia) con 31 años de registros mensuales (1993–2023), y está diseñado para su adaptación directa a cualquier estación de la red hidrométrica del IDEAM.

El marco integra dos módulos analíticos secuenciales:

- **Módulo 1 — Detección de no estacionariedad:** batería de seis pruebas formales de estacionariedad, análisis wavelet de Morlet multi-escala con simulación Monte Carlo, y descomposición de tendencias por época hidrológica.
- **Módulo 2 — Análisis de frecuencias:** ajuste estacionario de seis distribuciones de probabilidad con selección por AICc, y modelado no estacionario GAMLSS con validación por residuos de Cox-Snell e intervalos de confianza bootstrap.

El hallazgo central para el caso de estudio del Guachicono es una **paradoja estacionario/no estacionario** en *T*r ≈ 5–7 años: el modelo GEV clásico subestima los caudales de corto período de retorno en hasta +37,7 % y sobreestima los de largo período en hasta −70,5 %, respecto al modelo GAMLSS no estacionario para el año 2023. Este resultado tiene implicaciones directas para el dimensionamiento de infraestructura hidráulica en la región.

---

## Estructura del repositorio

```
CFE-Colombia/
│
├── CFE_Colombia_v1_0_0.R       # Script principal de análisis (1.284 líneas)
│
├── data/
│   └── Qmax.xlsx               # Datos de entrada: caudales máximos mensuales
│                               # (13 columnas: Año + 12 meses, m³/s)
│
├── outputs/                    # Carpeta generada automáticamente en la primera ejecución
│   ├── Figura1_Panel_Exploratorio.png
│   ├── Figura2_Test_Sneyers.png
│   ├── Figura3_Tendencia_Epoca.png
│   ├── Figura4_Extremos_Ventanas.png
│   ├── Figura5_Wavelet.png
│   ├── Figura6_Ajuste_Estacionario.png
│   ├── Figura7_Validacion_NoEstacionario.png
│   ├── Figura8_Cuantiles_NoEstacionarios.png
│   └── Informe_Final.txt
│
├── README.md                   # Esta documentación
└── LICENSE                     # CC BY 4.0
```

---

## Formato de los datos de entrada

El archivo `Qmax.xlsx` debe contener los **caudales máximos mensuales** en formato tabular ancho:

| Año  | Enero | Febrero | Marzo | … | Diciembre |
|------|-------|---------|-------|---|-----------|
| 1993 | 215,3 | 180,1   | 310,5 | … | 245,8     |
| 1994 | …     | …       | …     | … | …         |

- **Hoja:** 1 (primera hoja del libro)
- **Columnas:** 13 — `Año` (entero) seguido de los 12 meses en español (`Enero` a `Diciembre`)
- **Unidades:** m³/s
- **Datos faltantes:** celdas vacías; el script filtra los valores `NA` automáticamente

---

## Instalación

### Requisitos previos

- R ≥ 4.3.1
- RStudio (recomendado) o cualquier entorno compatible con R

### Paquetes requeridos

El script instala automáticamente los paquetes faltantes en la primera ejecución. La lista completa de dependencias es:

| Paquete | Versión | Función principal |
|---------|---------|-------------------|
| `readxl` | 1.4.5 | Lectura de datos Excel |
| `dplyr` | 1.1.4 | Manipulación de datos |
| `lubridate` | — | Manejo de fechas |
| `ggplot2` | 4.0.0 | Visualización |
| `tidyr` | — | Reestructuración de datos |
| `gridExtra` | 2.3 | Paneles de figuras múltiples |
| `gamlss` | 5.5-0 | Modelos GAMLSS no estacionarios |
| `gamlss.dist` | 6.1-1 | Familias de distribuciones GAMLSS |
| `gamlss.add` | — | Complementos GAMLSS |
| `ismev` | — | Ajuste GEV/GPD |
| `extRemes` | 2.2-1 | Bootstrap paramétrico GEV |
| `fitdistrplus` | 1.2-4 | Ajuste MLE (LNO, NOR) |
| `MASS` | — | Métodos estadísticos |
| `scales` | — | Escalas gráficas |
| `trend` | 1.1.6 | Pendiente de Sen, test de Pettitt |
| `Kendall` | 2.2.2 | Prueba de Mann-Kendall |
| `moments` | 0.14.1 | Momentos estadísticos |
| `evd` | 2.3-7.1 | Funciones Gumbel (`evd::dgumbel`) |
| `WaveletComp` | 1.2 | Análisis wavelet de Morlet |
| `zoo` | — | Estadísticas en ventanas móviles |
| `tseries` | 0.10-58 | Prueba ADF (raíz unitaria) |
| `lmtest` | 0.9-40 | Prueba de Breusch-Pagan (White) |
| `janitor` | — | Limpieza de datos |
| `grid` | — | Gráficos base (incluido en R) |

---

## Inicio rápido

### Paso 1 — Clonar el repositorio

```bash
git clone https://github.com/MauricioVictoriaN/CFE-Colombia.git
cd CFE-Colombia
```

### Paso 2 — Configurar las rutas

Abrir `CFE_Colombia_v1_0_0.R` y editar las dos variables de ruta en las líneas 257–258:

```r
ruta_archivo      <- "D:/R/Version_Final_CFE_Colombia/Qmax.xlsx"  # ← ruta al archivo de datos
directorio_salida <- "D:/R/Version_Final_CFE_Colombia/outputs"    # ← ruta para los resultados
```

Reemplazar ambas rutas con las ubicaciones reales en su sistema. La carpeta `outputs/` se crea automáticamente si no existe.

### Paso 3 — Ejecutar el script

En RStudio, abrir `CFE_Colombia_v1_0_0.R` y presionar **Ctrl + Alt + R** (Ejecutar todo), o desde la consola:

```r
source("CFE_Colombia_v1_0_0.R")
```

Una ejecución exitosa imprime en consola una lista de verificación de los archivos generados y guarda 9 archivos en `outputs/`.

---

## Metodología

### Módulo 1 — Detección de no estacionariedad

| Paso | Método | Criterio de decisión |
|------|--------|---------------------|
| 1.1 | Panel exploratorio (LOESS, boxplot estacional, evolución decadal) | Visual |
| 1.2 | Prueba de Ljung-Box de independencia serial (*k* = 6) | *p* > 0,05 para continuar |
| 1.3 | Mann-Kendall + pendiente de Sen | ≥ 2 / 4 pruebas formales significativas → no estacionaria |
| 1.4 | Prueba de cambio estructural de Pettitt | |
| 1.5 | Prueba de heterocedasticidad de White (Breusch-Pagan) | |
| 1.6 | Dickey-Fuller Aumentado (raíz unitaria) | |
| 1.7 | Sneyers secuencial (MK progresivo/retroactivo) | Complementaria |
| 1.8 | MK en ventanas móviles (*w* = 10 años) | Complementaria |
| 1.9 | Frecuencia de extremos en ventanas móviles | Visual |
| 1.10 | Espectrograma wavelet de Morlet (*n*_sim = 1.000 Monte Carlo) | Significancia al 5 % |

El criterio de enrutamiento sigue a Villarini et al. (2009): la serie se clasifica como **no estacionaria** si al menos 2 de las 4 pruebas con *p*-valor formal (MK, Pettitt, White, ADF) resultan significativas.

### Módulo 2 — Análisis de frecuencias

#### Rama estacionaria

Se ajustan seis distribuciones por Máxima Verosimilitud (MLE):

| ID | Distribución | Parámetros | Notas de implementación |
|----|-------------|------------|------------------------|
| GEV | Valor Extremo Generalizado | ξ, σ, κ | `extRemes` |
| GUM | Gumbel (VE Tipo I) | μ, σ | `evd::dgumbel` (evita conflictos de namespace) |
| LNO | Log-Normal | μ, σ | `fitdistrplus` |
| PE3 | Pearson III | α, β, ξ | Optimizador L-BFGS-B con límites explícitos |
| LP3 | Log-Pearson III | — | Inicialización por momentos logarítmicos |
| NOR | Normal | μ, σ | `fitdistrplus` |

Selección del mejor modelo por **AICc** (Burnham & Anderson 2002); IC 90 % mediante **bootstrap paramétrico** (*B* = 1.000) sobre GEV.

#### Rama no estacionaria (GAMLSS)

Se evalúan seis modelos:

| Modelo | Familia | μ(t) | σ(t) |
|--------|---------|------|------|
| M0 | Normal | constante | constante |
| M1 | Gumbel | lineal en *t* | constante |
| M2 | Gumbel | spline P | constante |
| M3 | Gumbel | spline P | lineal en *t* |
| **M4** | **Log-Normal** | **spline P** | **spline P** |
| M5 | Gamma | spline P | constante |

Selección por **AIC**. Validación mediante **residuos de Cox-Snell** (prueba KS + coeficiente de Filliben). IC 95 % mediante **bootstrap de residuos** (*B* = 500).

> **Nota de implementación:** La variable temporal se centra (`Año_centrado = Año − media(Año)`) para mejorar la estabilidad numérica del ajuste de los splines P. El bootstrap GAMLSS usa `gamlss()` directo en cada réplica (evitando `update()` que presenta problemas de resolución de entorno con splines).

---

## Parámetros clave

Todos los parámetros configurables están centralizados en la **Sección 0.3** del script (líneas 232–254):

```r
SEMILLA_GLOBAL   <- 2024L      # Semilla aleatoria global (reproducibilidad)
N_SIM_WAVELET    <- 1000L      # Simulaciones Monte Carlo para wavelet
B_BOOT_EST       <- 1000L      # Réplicas bootstrap — IC GEV estacionario (90 %)
B_BOOT_GAMLSS    <- 500L       # Réplicas bootstrap — IC GAMLSS (95 %)
ALPHA_BOOT       <- 0.10       # 1 − nivel de confianza para IC estacionario
DELTA_AIC_UMBRAL <- 2.0        # Umbral ΔAICc para competencia de modelos
METODO_TRAZADO   <- "cunnane"  # Método de posición de trazado
nivel_significancia <- 0.05    # α para todas las pruebas de hipótesis

# Períodos de retorno de diseño (años)
Tr_diseno <- c(2, 2.33, 5, 10, 15, 20, 25, 50, 100, 200, 500)

# Configuración de figuras
CFG_GRAF <- list(dpi = 300L, ancho_s = 10, ancho_d = 14, alto_p = 12)
```

---

## Archivos de salida

Todos los archivos se guardan en `directorio_salida` (por defecto: `outputs/`):

| Archivo | Descripción |
|---------|-------------|
| `Figura1_Panel_Exploratorio.png` | Panel 4-en-1: serie temporal completa, patrón estacional, evolución decadal, tendencia media anual |
| `Figura2_Test_Sneyers.png` | Test de Mann-Kendall Sneyers secuencial (curvas progresiva y retroactiva) |
| `Figura3_Tendencia_Epoca.png` | Tendencia de caudal medio por época hidrológica (Dic-Feb, Mar-May, Jun-Ago, Sep-Nov) |
| `Figura4_Extremos_Ventanas.png` | Frecuencia de excedencia de umbrales P90/P95/P99 en ventanas móviles de 10 años |
| `Figura5_Wavelet.png` | Espectrograma wavelet de Morlet con contornos de significancia Monte Carlo |
| `Figura6_Ajuste_Estacionario.png` | Ajuste de distribuciones estacionarias (comparación FDA, posición de Cunnane) |
| `Figura7_Validacion_NoEstacionario.png` | Q-Q plot de Cox-Snell del modelo GAMLSS ganador |
| `Figura8_Cuantiles_NoEstacionarios.png` | Evolución temporal de cuantiles no estacionarios (5 paneles × 5 períodos de retorno) con IC bootstrap |
| `Informe_Final.txt` | Informe numérico completo con todos los resultados y `sessionInfo()` |

---

## Resultados del caso de estudio — Río Guachicono (1993–2023)

| Indicador | Valor |
|-----------|-------|
| Longitud del registro | 31 años (372 registros mensuales) |
| Tendencia MK | τ = 0,316; pendiente de Sen = 2,48 m³/s/año; *p* = 0,013 |
| Cambio estructural (Pettitt) | 2010; Δ = 37,5 %; *p* = 0,020 |
| Mejor modelo estacionario | Pearson III (AICc = 374,68) |
| Mejor modelo GAMLSS | M4: Log-Normal con μ(t) y σ(t) splines P (AIC = 360,78) |
| Validación GAMLSS | KS *p* = 0,744; Filliben *r* = 0,976; edf_μ = 4,59; df.fit = 9,80 |
| Caudal de diseño *T*r = 100 años (2023) | 542,9 m³/s [IC 95 %: 455,2–2.435,2] |
| GEV estacionario *T*r = 100 años | 1.160,8 m³/s (−53,2 % vs. GAMLSS) |
| **Paradoja del cruce** | **Inversión en *T*r ≈ 5–7 años** |

### Resumen de la paradoja estacionario/no estacionario

| *T*r (años) | GEV estacionario (m³/s) | GAMLSS 2023 (m³/s) | Diferencia |
|-------------|------------------------|---------------------|------------|
| 2 | 324,0 | 446,1 | +37,7 % |
| 5 | 450,6 | 478,9 | +6,3 % |
| **~5–7** | **— punto de cruce —** | | |
| 10 | 563,9 | 497,1 | −11,8 % |
| 50 | 934,3 | 530,5 | −43,2 % |
| 100 | 1.160,8 | 542,9 | −53,2 % |
| 500 | 1.927,6 | 568,7 | −70,5 % |

---

## Adaptación a otras estaciones IDEAM

Para aplicar CFE-Colombia a una estación diferente, solo se requieren **dos modificaciones**:

1. Reemplazar `Qmax.xlsx` con los datos de caudal máximo mensual de la estación objetivo en el mismo formato de 13 columnas.
2. Actualizar `ruta_archivo` y `directorio_salida` con las rutas correspondientes al nuevo archivo.

Todos los parámetros de análisis (iteraciones bootstrap, nivel de significancia, períodos de retorno, método de posición de trazado) están centralizados en la **Sección 0.3** del script y pueden ajustarse de forma independiente sin modificar el resto del código.

---

## Reproducibilidad

El análisis es completamente reproducible:

- La semilla aleatoria global `SEMILLA_GLOBAL <- 2024L` está fijada al inicio del script mediante `set.seed()`.
- El informe de salida `Informe_Final.txt` incluye automáticamente el resultado de `sessionInfo()`, que registra las versiones exactas de R y todos los paquetes utilizados en cada ejecución.

---

## Cómo citar

Si utiliza este marco en su investigación, por favor cite el preprint asociado:

> Victoria Niño, M.J. (2026). *Análisis de frecuencias no estacionario de caudales máximos en el río Guachicono (Andes colombianos): detección de tendencias, modelado GAMLSS y cuantiles de diseño*. EngrXiv. https://doi.org/10.XXXX/engrxiv.XXXXXXXX

**BibTeX:**

```bibtex
@techreport{victorianino2026cfe,
  author       = {Victoria Ni{\~n}o, Mauricio Javier},
  title        = {An{\'a}lisis de frecuencias no estacionario de caudales m{\'a}ximos
                  en el r{\'i}o {Guachicono} ({Andes} colombianos): detecci{\'o}n de
                  tendencias, modelado {GAMLSS} y cuantiles de dise{\~n}o},
  institution  = {EngrXiv},
  year         = {2026},
  doi          = {10.XXXX/engrxiv.XXXXXXXX},
  url          = {https://github.com/MauricioVictoriaN/CFE-Colombia}
}
```

---

## Referencias

El marco implementa métodos de las siguientes referencias clave:

- Burnham, K.P. & Anderson, D.R. (2002). *Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach*. Springer. https://doi.org/10.1007/b97636
- Hamed, K.H. & Rao, A.R. (1998). A modified Mann-Kendall trend test for autocorrelated data. *Journal of Hydrology*, 204(1–4), 182–196. https://doi.org/10.1016/S0022-1694(97)00125-X
- Khaliq, M.N., Ouarda, T.B.M.J., Ondo, J.-C., Gachon, P. & Bobée, B. (2006). Frequency analysis of a sequence of dependent and/or non-stationary hydro-meteorological observations: A review. *Journal of Hydrology*, 329(3–4), 534–552. https://doi.org/10.1016/j.jhydrol.2006.03.004
- Poveda, G. & Álvarez, D.M. (2012). El colapso de la hipótesis de estacionariedad por cambio y variabilidad climática: implicaciones para el diseño hidrológico en ingeniería. *Revista de Ingeniería*, núm. 36, pp. 65–76. Universidad de los Andes. https://www.redalyc.org/articulo.oa?id=121025826012
- Rigby, R.A. & Stasinopoulos, D.M. (2005). Generalized additive models for location, scale and shape. *Applied Statistics*, 54(3), 507–554. https://doi.org/10.1111/j.1467-9876.2005.00510.x
- Salas, J.D. & Obeysekera, J. (2014). Revisiting the concepts of return period and risk for nonstationary hydrologic extreme events. *Journal of Hydrologic Engineering*, 19(3), 554–568. https://doi.org/10.1061/(ASCE)HE.1943-5584.0000820
- Villarini, G., Serinaldi, F., Smith, J.A. & Krajewski, W.F. (2009). On the stationarity of annual flood peaks in the continental United States during the 20th century. *Water Resources Research*, 45, W08417. https://doi.org/10.1029/2008WR007645

---

## Licencia

Este trabajo está licenciado bajo una [Licencia Creative Commons Atribución 4.0 Internacional](https://creativecommons.org/licenses/by/4.0/deed.es).

Es libre de compartir y adaptar el material para cualquier propósito, siempre que se otorgue el crédito apropiado, se incluya un enlace a la licencia y se indiquen los cambios realizados.


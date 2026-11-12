# =============================================================================
# STAR POR TESELAS + AMENAZA DOMINANTE (VERSIÓN DEFINITIVA Y ROBUSTA)
# =============================================================================

# Cargar paquetes necesarios para procesamiento espacial, tabular y visualización
library(terra)    # Procesamiento de datos ráster y vectoriales
library(sf)       # Manejo de capas vectoriales (shapefiles)
library(dplyr)    # Manipulación y resumen de tablas de datos
library(readr)    # Lectura rápida de archivos CSV
library(ggplot2)  # Creación de gráficos y mapas temáticos

# Verificar si el paquete 'exactextractr' está instalado; si no, instalarlo y cargarlo
if(!require(exactextractr)){
  install.packages("exactextractr")
  library(exactextractr)
}

# =============================================================================
# 1. RUTAS Y CONFIGURACIÓN
# =============================================================================

# Definir carpeta de entrada con los mapas ráster de Área de Hábitat (AOH)
dir_aoh      <- "C:/Users/walter.garcia/Documents/STAR/AOH_Generados_MM"
# Definir ruta de la tabla CSV con metadatos de las especies
csv_file     <- "C:/Users/walter.garcia/Documents/STAR/especies_mm2.csv"
# Definir ruta del shapefile con las teselas (unidades espaciales de análisis)
teselas_file <- "C:/Users/walter.garcia/Downloads/Teselas_Paisajes Focales Núcleos-20260409T192519Z-3-001/Teselas_Paisajes Focales Núcleos/Teselas_Nucleo_Magdalena_Medio/NMM_PF_Buffer_10km.shp"
# Definir carpeta que contiene los rásters de amenazas/presiones antrópicas
dir_pressure <- "C:/Users/walter.garcia/Documents/STAR/Presiones"
# Definir carpeta de salida donde se guardarán los mapas y shapefiles resultantes
output_dir   <- "C:/Users/walter.garcia/Documents/STAR/Resultados_STAR_Teselas_MM"

# Crear la carpeta de salida (y sus subcarpetas) si aún no existe
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# 2. CARGAR METADATA Y PESOS IUCN
# =============================================================================

# Leer el CSV de especies delimitado por punto y coma
meta <- read_delim(csv_file, delim = ";", show_col_types = FALSE)

# Definir tabla de conversión de categorías UICN a pesos numéricos STAR
iucn_weights <- c(
  "LC" = 0,    # Preocupación Menor
  "NT" = 100,  # Casi Amenazada
  "VU" = 200,  # Vulnerable
  "EN" = 300,  # En Peligro
  "CR" = 400   # En Peligro Crítico
)

# Asignar el peso correspondiente a cada fila según su categoría UICN
meta$weight <- iucn_weights[meta$iucn_status]

# Estandarizar nombres de especies: quitar espacios extra y pasar a minúsculas
meta$Species_clean <- trimws(tolower(meta$Species))

# =============================================================================
# 3. ARCHIVOS AOH
# =============================================================================

# Obtener la lista completa de rutas de archivos que terminan en "_AOH.tif"
aoh_files <- list.files(dir_aoh, pattern = "_AOH\\.tif$", full.names = TRUE)

# Detener el script si no se encuentra ningún archivo ráster
if(length(aoh_files) == 0){
  stop("No se encontraron rasters AOH en la carpeta especificada.")
}

# Cargar el primer ráster como plantilla de referencia para sistemas de coordenadas (CRS)
ref_raster <- rast(aoh_files[1])

# =============================================================================
# 4. CARGAR TESELAS Y REPROYECTAR
# =============================================================================

# Cargar el shapefile de polígonos (teselas)
teselas <- st_read(teselas_file, quiet = TRUE)
# Reproyectar las teselas al mismo CRS del ráster de referencia
teselas <- st_transform(teselas, crs(ref_raster))
# Crear un identificador único secuencial para cada tesela
teselas$ID <- 1:nrow(teselas)

# Convertir el objeto 'sf' a 'SpatVector' de terra para validaciones espaciales rápidas
teselas_vect <- vect(teselas)

# =============================================================================
# 5. CARGAR CAPAS DE PRESIÓN
# =============================================================================

# Listar todos los archivos .tif dentro de la carpeta de presiones
pressure_files <- list.files(dir_pressure, pattern = "\\.tif$", full.names = TRUE)

# Detener ejecución si no hay capas de presión
if(length(pressure_files) == 0){
  stop("No se encontraron capas de presión.")
}

# Inicializar lista vacía para almacenar los rásters de presión
pressure_rasters <- list()

# Recorrer cada capa de presión, reproyectarla si es necesario y guardarla en la lista
for(f in pressure_files){
  # Extraer el nombre del archivo sin extensión y convertirlo en un nombre válido en R
  nombre <- tools::file_path_sans_ext(basename(f))
  nombre <- make.names(nombre)
  
  # Cargar el ráster de presión
  r <- rast(f)
  # Reproyectar al CRS de referencia si no coinciden
  if(!same.crs(r, ref_raster)){
    r <- project(r, ref_raster)
  }
  # Guardar en la lista usando el nombre limpio como llave
  pressure_rasters[[nombre]] <- r
}

# Definir función para escalar los valores de un ráster entre 0 y 1 (normalización Min-Max)
normalize_raster <- function(r){
  # Extraer valores válidos (no NA) y finitos
  vals <- values(r, na.rm = TRUE)
  vals <- vals[is.finite(vals)]
  # Si el ráster está vacío o es constante, devolver un ráster de ceros
  if(length(vals) == 0) return(r * 0)
  rmin <- min(vals)
  rmax <- max(vals)
  if(rmax == rmin) return(r * 0)
  # Aplicar fórmula de normalización: (x - min) / (max - min)
  return((r - rmin) / (rmax - rmin))
}

# Aplicar la normalización 0-1 a todas las capas de presión cargadas
pressure_norm <- lapply(pressure_rasters, normalize_raster)

# =============================================================================
# 6. INICIALIZAR CAMPOS DE SALIDA EN TESELAS
# =============================================================================

# Crear columnas en cero para ir acumulando el STAR total y el STAR por cada tipo de amenaza
teselas$STAR_TOTAL <- 0
teselas$STAR_FRAG  <- 0
teselas$STAR_LU    <- 0
teselas$STAR_POP   <- 0
teselas$STAR_VIAS  <- 0

# =============================================================================
# 7. BUCLE PRINCIPAL DE PROCESAMIENTO POR ESPECIE
# =============================================================================

cat("\n============================================================\n")
cat("PROCESANDO ESPECIES Y CALCULANDO STAR POR TESELA\n")
cat("============================================================\n")

# Iterar sobre cada archivo AOH de especie
for(i in seq_along(aoh_files)){
  
  aoh_path <- aoh_files[i]
  nombre_archivo <- basename(aoh_path)
  
  # --- Extracción con Expresión Regular (Género_especie) ---
  # Buscar un patrón de texto tipo "Genero_especie" dentro del nombre del archivo
  coincidencia <- regmatches(
    nombre_archivo, 
    regexpr("[A-Z][a-z]+_[a-z]+", nombre_archivo)
  )
  
  # Si el archivo no tiene un nombre científico reconocible, saltar al siguiente
  if(length(coincidencia) == 0){
    next
  }
  
  # Crear versiones del nombre con espacio ("Panthera onca") y con guion bajo ("Panthera_onca")
  spname_espacio <- gsub("_", " ", coincidencia)
  spname_snake   <- coincidencia
  
  # --- Buscar coincidencia en metadata ---
  # Filtrar el CSV buscando coincidencias en cualquiera de las columnas de nombre/archivo
  row_meta <- meta %>% 
    filter(
      Species_clean == trimws(tolower(spname_espacio)) |
        tolower(trimws(Species_2)) == tolower(spname_snake) |
        grepl(tolower(spname_snake), tolower(species_file), fixed = TRUE)
    )
  
  # Si la especie no aparece en el CSV, reportar en consola y omitir
  if(nrow(row_meta) == 0){
    cat(sprintf("[%d/%d] Especie: '%s' (de %s) -> Sin metadata en CSV\n", i, length(aoh_files), spname_espacio, nombre_archivo))
    next
  }
  
  # Obtener el peso UICN de la especie e imprimir avance en consola
  peso <- row_meta$weight[1]
  cat(sprintf("[%d/%d] Especie encontrada: %s (%s, Peso: %d)\n", i, length(aoh_files), row_meta$Species[1], row_meta$iucn_status[1], peso))
  
  # --- Cargar AOH de la especie ---
  r <- rast(aoh_path)
  # Convertir a NA cualquier celda que no sea 1 (presencia de hábitat)
  r[r != 1] <- NA
  
  # Calcular el número total de celdas de hábitat global de la especie
  global_cells <- global(!is.na(r), "sum", na.rm = TRUE)[1, 1]
  if(is.na(global_cells) || global_cells == 0){
    cat("  -> Sin presencia global (celdas vacías)\n")
    next
  }
  
  # Verificar rápidamente si la extensión del AOH intersecta con el área de las teselas
  if(!is.related(r, teselas_vect, "intersects")){
    cat("  -> Sin solape geográfico con el área de teselas\n")
    next
  }
  
  # --- Extracción en teselas ---
  # Sumar la fracción de celdas de hábitat dentro de cada polígono de manera vectorizada y rápida
  cov_list <- tryCatch({
    exact_extract(r, teselas, fun = "sum", progress = FALSE)
  }, error = function(e) NULL)
  
  if(is.null(cov_list)) next
  cov_list[is.na(cov_list)] <- 0
  
  # Si ninguna tesela contiene hábitat de esta especie, pasar a la siguiente
  if(sum(cov_list) == 0){
    cat("  -> Sin presencia dentro de las teselas de estudio\n")
    next
  }
  
  # Calcular proporción de hábitat por tesela (START) y multiplicarla por el peso UICN (STAR)
  start_vec <- cov_list / global_cells
  star_vec  <- start_vec * peso
  
  # Sumar la contribución STAR de esta especie al total acumulado de cada tesela
  teselas$STAR_TOTAL <- teselas$STAR_TOTAL + star_vec
  
  # --- Extraer Presiones re-alineando dinámicamente si difieren las geometrías ---
  for(p in names(pressure_norm)){
    rp <- pressure_norm[[p]]
    
    # Si la capa de presión no coincide en resolución o extensión con el AOH, remuestrear
    if(!same.crs(rp, r) || !compareGeom(rp, r, stopOnError = FALSE)){
      rp <- resample(rp, r, method = "bilinear")
    }
    
    # Calcular la presión promedio en cada tesela ponderada únicamente por las celdas donde está la especie
    p_mean <- tryCatch({
      exact_extract(rp, teselas, fun = "mean", weights = r, progress = FALSE)
    }, error = function(e) NULL)
    
    if(is.null(p_mean)) next
    p_mean[is.na(p_mean)] <- 0
    
    # Ponderar el STAR de la especie en cada tesela por la intensidad de esta presión específica
    star_p_vec <- star_vec * p_mean
    
    # Acumular el resultado en la columna correspondiente según el nombre de la capa de presión
    if(p == "frag2022")      teselas$STAR_FRAG <- teselas$STAR_FRAG + star_p_vec
    else if(p == "LU12022")  teselas$STAR_LU   <- teselas$STAR_LU   + star_p_vec
    else if(p == "Pop2022")  teselas$STAR_POP  <- teselas$STAR_POP  + star_p_vec
    else if(p == "Vias2022") teselas$STAR_VIAS <- teselas$STAR_VIAS + star_p_vec
  }
}

# =============================================================================
# 8. DETERMINAR AMENAZA DOMINANTE POR TESELA
# =============================================================================

cat("\nCalculando amenaza dominante por tesela...\n")

# Evaluar fila por fila (tesela por tesela) cuál de las 4 presiones aporta mayor valor STAR
teselas <- teselas %>%
  rowwise() %>%
  mutate(
    # Sumar los 4 componentes de amenaza
    Suma_Amenazas = sum(c(STAR_FRAG, STAR_LU, STAR_POP, STAR_VIAS), na.rm = TRUE),
    # Asignar etiqueta "Sin Amenaza" si la suma es 0; de lo contrario, elegir el nombre del valor máximo
    Amenaza_Dominante = if_else(
      Suma_Amenazas == 0,
      "Sin Amenaza / Presencia",
      c("Fragmentacion", "Uso_Suelo", "Poblacion", "Vias")[
        which.max(c(STAR_FRAG, STAR_LU, STAR_POP, STAR_VIAS))
      ]
    )
  ) %>%
  ungroup()

# =============================================================================
# 9. MAPA Y EXPORTACIÓN
# =============================================================================

# Definir paleta de colores hexadecimales para cada categoría de amenaza
colores_amenazas <- c(
  "Fragmentacion"           = "#d73027",
  "Uso_Suelo"               = "#fc8d59",
  "Poblacion"               = "#4575b4",
  "Vias"                    = "#1a9850",
  "Sin Amenaza / Presencia" = "#cccccc"
)

# Construir el mapa temático de teselas usando ggplot2
mapa <- ggplot(teselas) +
  geom_sf(aes(fill = Amenaza_Dominante), color = "white", linewidth = 0.05) +
  scale_fill_manual(values = colores_amenazas) +
  theme_minimal() +
  labs(
    title = "Amenaza Dominante por Tesela (STAR)",
    fill  = "Categoría"
  )

# Exportar el objeto espacial actualizado como Shapefile
st_write(
  teselas,
  file.path(output_dir, "STAR_Teselas_Amenaza_Dominante.shp"),
  delete_layer = TRUE
)

# Guardar el mapa generado como imagen PNG de alta resolución (300 dpi)
ggsave(
  file.path(output_dir, "Mapa_Amenaza_Dominante.png"),
  mapa,
  width = 10,
  height = 8,
  dpi = 300
)

cat("\n¡Proceso finalizado exitosamente!\nResultados en:", output_dir, "\n")

# --- Bloque de exploración y estadísticas descriptivas en consola ---
print(mapa)                                  # Mostrar el mapa en el visor de R
summary(teselas$STAR_POP)                    # Resumen estadístico de STAR por población
max(teselas$STAR_POP, na.rm = TRUE)          # Valor máximo de STAR por población
sum(teselas$STAR_POP > 0, na.rm = TRUE)      # Cantidad de teselas afectadas por población

summary(teselas$STAR_FRAG)                   # Resumen estadístico de STAR por fragmentación
summary(teselas$STAR_LU)                     # Resumen estadístico de STAR por uso de suelo
summary(teselas$STAR_VIAS)                   # Resumen estadístico de STAR por vías

table(teselas$Amenaza_Dominante)             # Conteo de teselas por amenaza dominante

# Calcular el porcentaje de teselas correspondiente a cada amenaza dominante
round(
  100 * prop.table(
    table(teselas$Amenaza_Dominante)
  ),
  2
)


# =============================================================================
# GENERAR TABLA DE PRESIONES PROMEDIO POR ESPECIE
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(readr)

# 1. Rutas
dir_aoh      <- "C:/Users/walter.garcia/Documents/STAR/AOH_Generados_MM"
csv_file     <- "C:/Users/walter.garcia/Documents/STAR/especies_mm2.csv"
dir_pressure <- "C:/Users/walter.garcia/Documents/STAR/Presiones"
output_dir   <- "C:/Users/walter.garcia/Documents/STAR/Resultados_MM_presiones"

# Crear directorio de salida para las tablas de especies si no existe
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. Cargar Metadata y Rasters
meta <- read_delim(csv_file, delim = ";", show_col_types = FALSE)
aoh_files <- list.files(dir_aoh, pattern = "_AOH\\.tif$", full.names = TRUE)
ref_raster <- rast(aoh_files[1])

# 3. Cargar y Normalizar Capas de Presión
pressure_files <- list.files(dir_pressure, pattern = "\\.tif$", full.names = TRUE)
pressure_norm <- list()

# Cargar, reproyectar y normalizar (0-1) cada ráster de presión
for(f in pressure_files){
  nombre <- tools::file_path_sans_ext(basename(f))
  nombre <- make.names(nombre)
  r <- rast(f)
  if(!same.crs(r, ref_raster)) r <- project(r, ref_raster)
  
  # Normalizar valores entre 0 y 1
  vals <- values(r, na.rm = TRUE)
  vals <- vals[is.finite(vals)]
  rmin <- min(vals); rmax <- max(vals)
  r_norm <- (r - rmin) / (rmax - rmin)
  
  pressure_norm[[nombre]] <- r_norm
}

# 4. Procesar Especie por Especie
lista_resultados <- list()

# Recorrer cada especie para calcular la presión promedio dentro de todo su hábitat
for(i in seq_along(aoh_files)){
  
  aoh_path <- aoh_files[i]
  nombre_archivo <- basename(aoh_path)
  
  # Extraer Nombre Científico desde el nombre del archivo
  coincidencia <- regmatches(nombre_archivo, regexpr("[A-Z][a-z]+_[a-z]+", nombre_archivo))
  if(length(coincidencia) == 0) next
  
  spname_espacio <- gsub("_", " ", coincidencia)
  spname_snake   <- coincidencia
  
  # Buscar Categoría IUCN en el CSV (si no existe asigna "NE/Sin Meta")
  row_meta <- meta %>% 
    filter(
      tolower(trimws(Species)) == tolower(spname_espacio) |
        tolower(trimws(Species_2)) == tolower(spname_snake)
    )
  
  cat_iucn <- if(nrow(row_meta) > 0) as.character(row_meta$iucn_status[1]) else "NE/Sin Meta"
  
  # Cargar ráster AOH y dejar solo píxeles con valor 1
  r <- rast(aoh_path)
  r[r != 1] <- NA
  
  # Omitir la especie si su mapa no tiene celdas de presencia
  if(global(!is.na(r), "sum", na.rm=TRUE)[1,1] == 0) next
  
  # Crear fila inicial de resultados con el nombre y categoría de la especie
  res_row <- data.frame(
    ESPECIE = spname_espacio,
    CATEGORIA = cat_iucn,
    stringsAsFactors = FALSE
  )
  
  # Extraer Promedios de Presión sobre el hábitat de la especie
  for(p in names(pressure_norm)){
    rp <- pressure_norm[[p]]
    
    # Alinear geometría de la capa de presión con la del AOH si difieren
    if(!compareGeom(rp, r, stopOnError = FALSE)){
      rp <- resample(rp, r, method = "bilinear")
    }
    
    # Extraer los valores de presión únicamente en los píxeles donde la especie está presente
    vals_p <- values(rp)[!is.na(values(r))]
    vals_p <- vals_p[!is.na(vals_p) & is.finite(vals_p)]
    
    # Calcular la media redondeada a 3 decimales (o 0 si no hay datos)
    val_mean <- if(length(vals_p) > 0) round(mean(vals_p), 3) else 0
    
    # Renombrar códigos de archivo a nombres de columna legibles
    col_name <- case_when(
      p == "frag2022" ~ "FRAGMENTACION",
      p == "LU12022"  ~ "USO_DEL_SUELO",
      p == "Pop2022"  ~ "POBLACION",
      p == "Vias2022" ~ "VIAS",
      TRUE ~ p
    )
    
    # Agregar el promedio calculado como nueva columna en la fila de la especie
    res_row[[col_name]] <- val_mean
  }
  
  # Añadir la fila completa de la especie a la lista general
  lista_resultados[[length(lista_resultados) + 1]] <- res_row
}

# 5. Unir y Exportar Tabla
# Consolidar todas las filas en una tabla única y eliminar posibles especies duplicadas
tabla_especies <- bind_rows(lista_resultados) %>% distinct(ESPECIE, .keep_all = TRUE)

# Mostrar las primeras 10 filas en consola para verificación
print(head(tabla_especies, 10))

# Exportar la tabla de presiones promedio absolutas (0-1) a CSV
write.csv(
  tabla_especies, 
  file.path(output_dir, "Tabla_Especies_Presiones_Completa.csv"), 
  row.names = FALSE
)

cat("\n¡Tabla generada con éxito en!:", file.path(output_dir, "Tabla_Especies_Presiones_Completa.csv"))


# =============================================================================
# 6. CONVERTIR PRESIONES DE ESPECIE A PORCENTAJES RELATIVOS
# =============================================================================

library(dplyr)

# Calcular qué porcentaje representa cada presión respecto a la presión total que sufre cada especie
tabla_porcentajes <- tabla_especies %>%
  rowwise() %>%
  mutate(
    # Sumar los 4 valores promedio de presión de la especie
    Suma_Total = sum(c(FRAGMENTACION, USO_DEL_SUELO, POBLACION, VIAS), na.rm = TRUE),
    
    # Calcular el porcentaje relativo de cada amenaza (redondeado a 1 decimal)
    Pct_Fragmentacion = round((FRAGMENTACION / Suma_Total) * 100, 1),
    Pct_Uso_Suelo     = round((USO_DEL_SUELO / Suma_Total) * 100, 1),
    Pct_Poblacion     = round((POBLACION / Suma_Total) * 100, 1),
    Pct_Vias          = round((VIAS / Suma_Total) * 100, 1)
  ) %>%
  ungroup() %>%
  # Conservar únicamente las columnas de identificación y los porcentajes calculados
  select(ESPECIE, CATEGORIA, Pct_Fragmentacion, Pct_Uso_Suelo, Pct_Poblacion, Pct_Vias)

# Mostrar las primeras filas de la tabla porcentual en consola
print(head(tabla_porcentajes))

# Exportar la tabla de porcentajes a un archivo CSV en el directorio de trabajo actual
write.csv(tabla_porcentajes, "Tabla_Especies_Presiones_Porcentaje.csv", row.names = FALSE)

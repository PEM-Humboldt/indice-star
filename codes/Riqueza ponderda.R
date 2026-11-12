# Cargar las librerías necesarias para análisis espacial y manipulación de datos
library(terra)  # Manejo y procesamiento de datos ráster y vectoriales
library(sf)     # Lectura, escritura y transformación de datos vectoriales (shapefiles)
library(dplyr)  # Manipulación, filtrado y resumen de tablas de datos
library(readr)  # Lectura rápida y eficiente de archivos CSV

# ===============================
# 1. RUTAS
# ===============================
# Definir la carpeta donde se encuentran los mapas ráster de Área de Hábitat (AOH)
dir_aoh <- "C:/Users/walter.garcia/Documents/STAR/AOH_Generados_3"

# Definir la ruta del archivo CSV con la información taxonómica y de amenaza
csv_file <- "C:/Users/walter.garcia/Documents/STAR/Especies_amenazadas.csv"

# Definir la ruta del shapefile con los polígonos o unidades de análisis (teselas)
aoi_file <- "C:/Users/walter.garcia/Documents/STAR/STAR_test_10_density.shp"

# ===============================
# 2. METADATA
# ===============================
# Leer el archivo CSV separado por punto y coma sin mostrar mensajes de tipo de columna
meta <- read_delim(
  csv_file,
  delim = ";",
  show_col_types = FALSE
)

# Extraer únicamente el nombre del archivo (sin la ruta completa) para poder cruzarlo después
meta$species_file <- basename(meta$species_file)

# ===============================
# 3. PESOS IUCN
# ===============================
# Crear un vector nombrado con los puntajes asignados a cada categoría de amenaza UICN
iucn_weights <- c(
  "LC" = 0,    # Preocupación Menor (se excluye del cálculo)
  "NT" = 100,  # Casi Amenazada
  "VU" = 200,  # Vulnerable
  "EN" = 300,  # En Peligro
  "CR" = 400   # En Peligro Crítico
)

# Asignar el peso correspondiente a cada especie según su estado UICN en la tabla
meta$weight <- iucn_weights[meta$iucn_status]

# ===============================
# 4. TESSELAS
# ===============================
# Cargar el shapefile de las teselas en memoria sin imprimir detalles en consola
teselas <- st_read(aoi_file, quiet = TRUE)

# ===============================
# 5. RASTERS
# ===============================
# Listar todos los archivos que terminen exactamente en "_AOH.tif" dentro de la carpeta
r_files <- list.files(
  dir_aoh,
  pattern = "_AOH\\.tif$",
  full.names = TRUE
)

# Verificar si se encontraron archivos; si la lista está vacía, detener la ejecución con error
if (length(r_files) == 0) {
  stop("❌ No hay AOHs")
}

# ===============================
# 6. RASTER REFERENCIA
# ===============================
# Cargar el primer ráster de la lista para usarlo como plantilla de proyección y resolución
ref <- rast(r_files[1])

# ===============================
# 7. CRS TESSELAS
# ===============================
# Reproyectar las teselas al mismo Sistema de Referencia de Coordenadas (CRS) del ráster base
teselas <- st_transform(
  teselas,
  crs(ref)
)

# Convertir el objeto vectorial de 'sf' a formato 'SpatVector' de 'terra' para mayor velocidad
teselas_vect <- vect(teselas)

# ===============================
# 8. PRECALCULAR ÁREA GLOBAL
# ===============================
# Imprimir mensaje de progreso en consola
cat("📦 Calculando áreas globales...\n")

# Crear una lista vacía para almacenar el conteo total de celdas de cada especie
global_stats <- list()

# Iniciar bucle para recorrer cada archivo ráster de especie
for (f in r_files) {
  
  # Cargar el ráster de la especie actual
  r <- rast(f)
  
  # ---------------------------
  # ALINEAR SI NECESARIO
  # ---------------------------
  # Comprobar si el ráster actual coincide en proyección, extensión y resolución con 'ref'
  if (!compareGeom(r, ref, stopOnError = FALSE)) {
    
    # Reproyectar el ráster a la cuadrícula de 'ref' usando vecino más cercano ("near")
    r <- project(
      r,
      ref,
      method = "near"
    )
    
    # Asegurar que los píxeles estén exactamente alineados celda por celda con 'ref'
    r <- resample(
      r,
      ref,
      method = "near"
    )
  }
  
  # ---------------------------
  # BINARIO
  # ---------------------------
  # Convertir el mapa a binario: mantener los píxeles con valor 1 y convertir el resto en NA
  r_bin <- ifel(r == 1, 1, NA)
  
  # ---------------------------
  # GLOBAL CELLS
  # ---------------------------
  # Sumar todos los píxeles con valor 1 para obtener el área total global de la especie
  global_cells <- global(
    r_bin,
    "sum",
    na.rm = TRUE
  )[1,1]
  
  # Guardar el resultado en la lista usando el nombre del archivo como identificador
  global_stats[[basename(f)]] <- global_cells
}

# ===============================
# 9. LOOP PRINCIPAL
# ===============================
# Pre-asignar una lista vacía con el tamaño exacto del número de teselas
resultados_lista <- vector(
  "list",
  length = nrow(teselas_vect)
)

# Bucle externo: recorrer cada polígono (tesela) uno por uno
for (i in seq_len(nrow(teselas_vect))) {
  
  # Extraer el polígono actual
  poly <- teselas_vect[i]
  
  # Asignar un ID numérico secuencial a la tesela
  poly_id <- i
  
  # Mostrar el avance por tesela en la consola
  cat(
    "\n🧩 Tesela",
    i,
    "/",
    nrow(teselas_vect),
    "\n"
  )
  
  # Crear una lista temporal para guardar los resultados de todas las especies en esta tesela
  resultados_poly <- list()
  
  # ---------------------------
  # LOOP ESPECIES
  # ---------------------------
  # Bucle interno: evaluar cada especie dentro de la tesela actual
  for (j in seq_along(r_files)) {
    
    # Obtener ruta y nombre base del archivo ráster de la especie
    r_file <- r_files[j]
    spname_file <- basename(r_file)
    
    # ---------------------------
    # METADATA
    # ---------------------------
    # Buscar la fila correspondiente a esta especie en la tabla de metadatos
    row_meta <- meta %>%
      filter(species_file == spname_file)
    
    # Si la especie no está en el CSV, saltar a la siguiente iteración
    if (nrow(row_meta) == 0) next
    
    # Extraer el peso UICN de la especie
    peso <- row_meta$weight[1]
    
    # Si el peso es nulo (NA) o cero (categoría LC), omitir esta especie
    if (is.na(peso) || peso == 0) next
    
    # ---------------------------
    # RASTER
    # ---------------------------
    # Cargar el ráster de la especie desde el disco
    r <- rast(r_file)
    
    # ---------------------------
    # ALINEAR SI NECESARIO
    # ---------------------------
    # Verificar alineación geométrica con el ráster de referencia
    if (!compareGeom(r, ref, stopOnError = FALSE)) {
      
      # Reproyectar si hay diferencias
      r <- project(
        r,
        ref,
        method = "near"
      )
      
      # Remuestrear para ajustar celdas
      r <- resample(
        r,
        ref,
        method = "near"
      )
    }
    
    # ---------------------------
    # BINARIO
    # ---------------------------
    # Binarizar el ráster (1 = presencia de hábitat, NA = ausencia)
    r_bin <- ifel(r == 1, 1, NA)
    
    # ---------------------------
    # GLOBAL
    # ---------------------------
    # Recuperar el total de celdas globales precalculado en la Sección 8
    global_cells <- global_stats[[spname_file]]
    
    # Si la especie no tiene celdas globales válidas, saltar a la siguiente
    if (is.na(global_cells) || global_cells == 0) next
    
    # ---------------------------
    # CROP + MASK
    # ---------------------------
    # Recortar el ráster al rectángulo delimitador (bounding box) de la tesela actual
    r_crop <- crop(r_bin, poly)
    
    # Enmascarar el recorte para conservar únicamente los píxeles dentro del contorno exacto del polígono
    r_mask <- mask(r_crop, poly)
    
    # ---------------------------
    # OVERLAP
    # ---------------------------
    # Contar cuántos píxeles de hábitat de la especie caen dentro de la tesela
    overlap_cells <- global(
      r_mask,
      "sum",
      na.rm = TRUE
    )[1,1]
    
    # Si el resultado es NA (sin superposición), convertirlo a 0
    overlap_cells <- ifelse(
      is.na(overlap_cells),
      0,
      overlap_cells
    )
    
    # ---------------------------
    # PRESENCIA
    # ---------------------------
    # Determinar presencia binaria: 1 si hay al menos un píxel dentro de la tesela, 0 si no
    presence <- ifelse(
      overlap_cells > 0,
      1,
      0
    )
    
    # ---------------------------
    # START
    # ---------------------------
    # Calcular la proporción del hábitat global de la especie que está dentro de esta tesela
    START <- overlap_cells / global_cells
    
    # ---------------------------
    # STAR
    # ---------------------------
    # Multiplicar la proporción de hábitat por el peso de amenaza UICN de la especie
    STAR_species <- START * peso
    
    # Almacenar todas las métricas calculadas en una fila de tabla (data.frame)
    resultados_poly[[j]] <- data.frame(
      GRID_ID = poly_id,
      species = row_meta$Species[1],
      file = spname_file,
      iucn_status = row_meta$iucn_status[1],
      weight = peso,
      global_cells = global_cells,
      overlap_cells = overlap_cells,
      presence = presence,
      START = START,
      STAR_species = STAR_species
    )
  }
  
  # Unir los resultados de todas las especies evaluadas en esta tesela en una sola tabla
  resultados_lista[[i]] <- bind_rows(resultados_poly)
}

# ===============================
# 10. UNIR RESULTADOS
# ===============================
# Consolidar las tablas de todas las teselas en un único data.frame maestro
resultados_final <- bind_rows(resultados_lista)

# ===============================
# 11. STAR TOTAL
# ===============================
# Agrupar por tesela (GRID_ID) y sumar los valores STAR individuales de todas sus especies
STAR_por_tesela <- resultados_final %>%
  group_by(GRID_ID) %>%
  summarise(
    STAR_total = sum(
      STAR_species,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ===============================
# 12. RIQUEZA SIMPLE
# ===============================
# Calcular el número total de especies amenazadas/casi amenazadas presentes por tesela
riqueza_simple <- resultados_final %>%
  group_by(GRID_ID) %>%
  summarise(
    richness = sum(
      presence,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ===============================
# 13. RIQUEZA PONDERADA
# ===============================
# Sumar los pesos UICN de las especies presentes en cada tesela (sin considerar proporción de área)
riqueza_ponderada <- resultados_final %>%
  group_by(GRID_ID) %>%
  summarise(
    richness_weighted = sum(
      weight * presence,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ===============================
# 14. FUSIONAR MÉTRICAS
# ===============================
# Unir las tres tablas resumen (STAR, riqueza simple y riqueza ponderada) usando GRID_ID
fusion <- STAR_por_tesela %>%
  left_join(
    riqueza_simple,
    by = "GRID_ID"
  ) %>%
  left_join(
    riqueza_ponderada,
    by = "GRID_ID"
  )

# ===============================
# 15. NORMALIZAR
# ===============================
# Escalar STAR_total entre 0 y 1 dividiendo por el valor máximo de toda la grilla
fusion$STAR_norm <- (
  fusion$STAR_total /
    max(fusion$STAR_total, na.rm = TRUE)
)

# Escalar la riqueza simple entre 0 y 1
fusion$richness_norm <- (
  fusion$richness /
    max(fusion$richness, na.rm = TRUE)
)

# Escalar la riqueza ponderada entre 0 y 1
fusion$richness_weighted_norm <- (
  fusion$richness_weighted /
    max(fusion$richness_weighted, na.rm = TRUE)
)

# ===============================
# 16. ÍNDICE HÍBRIDO
# ===============================
# Calcular el índice compuesto dando 70% de importancia a STAR y 30% a riqueza ponderada
fusion$STAR_RICH <- (
  fusion$STAR_norm * 0.7 +
    fusion$richness_weighted_norm * 0.3
)

# ===============================
# 17. UNIR AL SHAPE
# ===============================
# Crear la columna identificadora GRID_ID en el objeto espacial original
teselas$GRID_ID <- 1:nrow(teselas)

# Unir todas las métricas calculadas a la geometría espacial de las teselas
teselas_STAR <- teselas %>%
  left_join(
    fusion,
    by = "GRID_ID"
  )

# ===============================
# 18. EXPORTAR
# ===============================
# Guardar el mapa vectorial resultante en formato Shapefile (sobrescribiendo si ya existe)
st_write(
  teselas_STAR,
  "C:/Users/walter.garcia/Documents/STAR/STAR_RICH_teselas.shp",
  delete_dsn = TRUE,
  quiet = TRUE
)

# Guardar la tabla detallada (especie por tesela) en formato CSV sin nombres de fila
write.csv(
  resultados_final,
  "C:/Users/walter.garcia/Documents/STAR/STAR_RICH_detalle.csv",
  row.names = FALSE
)

# Imprimir mensaje de finalización exitosa en la consola
cat("\n✅ STAR + riqueza completado\n")



# mapa --------------------------------------------------------------------


library(sf)
library(ggplot2)

# Leer shapefile
teselas <- st_read("C:/Users/walter.garcia/Documents/STAR/STAR riqueza/STAR_RICH_teselas.shp")

# Ver nombres de campos
names(teselas)


library(ggplot2)
library(sf)

mapa <- ggplot(teselas) +
  geom_sf(
    aes(fill = rchnss_w),
    color = "black",
    linewidth = 0.1
  ) +
  scale_fill_viridis_c(
    name = "Valor"
  ) +
  theme_minimal()

print(mapa)

library(sf)
library(dplyr)
library(ggplot2)

teselas <- teselas %>%
  mutate(
    clase = cut(
      rchnss_w,
      breaks = c(0, 1600, 2200, 2600, 3600, 4600),
      include.lowest = TRUE,
      labels = c(
        "0 - 1600",
        "1600 - 2200",
        "2200 - 2600",
        "2600 - 3600",
        "3600 - 4600"
      )
    )
  )

mapa <- ggplot(teselas) +
  
  geom_sf(
    aes(fill = clase),
    color = "black",
    linewidth = 0.05
  ) +
  
  scale_fill_manual(
    values = c(
      "0 - 1600" = "#d80027",
      "1600 - 2200" = "#f4a582",
      "2200 - 2600" = "#f7f7f7",
      "2600 - 3600" = "#bdbdbd",
      "3600 - 4600" = "#4d4d4d"
    ),
    drop = FALSE
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "left",
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +
  
  labs(fill = NULL)

print(mapa)




library(grid)

mapa <- ggplot(teselas) +
  
  geom_sf(
    aes(fill = clase),
    color = "black",
    linewidth = 0.05
  ) +
  
  scale_fill_manual(
    values = c(
      "0 - 1600" = "#d80027",
      "1600 - 2200" = "#f4a582",
      "2200 - 2600" = "#f7f7f7",
      "2600 - 3600" = "#bdbdbd",
      "3600 - 4600" = "#4d4d4d"
    ),
    drop = FALSE,
    name = NULL
  ) +
  
  labs(
    title = "Riqueza ponderada por tesela"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(
      hjust = 0,
      size = 12,
      face = "plain"
    ),
    
    legend.position = "left",
    
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    
    panel.grid.major = element_line(
      colour = "grey85",
      linewidth = 0.3
    ),
    
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      colour = NA
    )
  )

print(mapa)


ggsave(
  "Mapa_riqueza.png",
  plot = mapa,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

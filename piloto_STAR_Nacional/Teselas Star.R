# =========================================
# STAR POR TESELAS (VERSIÓN OPTIMIZADA Y RÁPIDA)
# =========================================

library(terra)
library(sf)
library(dplyr)
library(readr)

# ===============================
# 1. RUTAS
# ===============================
dir_aoh  <- "C:/Users/walter.garcia/Documents/Piloto STAR Colombia/AOH"
csv_file <- "C:/Users/walter.garcia/Documents/STAR/Especies_amenazadas.csv"
aoi_file <- "C:/Users/walter.garcia/Documents/Piloto STAR Colombia/teselas_hexagonos_colombia_50km.shp"

# ===============================
# 2. METADATA
# ===============================
meta <- read_delim(csv_file, delim = ";", show_col_types = FALSE)
meta$species_file <- basename(meta$species_file)

iucn_weights <- c("LC" = 0, "NT" = 100, "VU" = 200, "EN" = 300, "CR" = 400)
meta$weight <- iucn_weights[meta$iucn_status]

# ===============================
# 3. TESELAS (AOI)
# ===============================
teselas <- st_read(aoi_file)
teselas$GRID_ID <- 1:nrow(teselas)


# Cargar el shapefile del límite de Colombia (o descargarlo si no lo tienes)
# Asegúrate de reprojectarlo al mismo CRS que las teselas
colombia_contorno <- st_read("C:/Users/walter.garcia/Documents/Indice de vacios/colombia_shapefile.shp") %>% 
  st_transform(crs = st_crs(teselas))

# Filtrar las teselas para quedarte SOLO con las que intersecan la tierra
teselas <- st_filter(teselas, colombia_contorno)

# Re-asignar IDs limpios
teselas$GRID_ID <- 1:nrow(teselas)
teselas_proj <- st_transform(teselas, crs = ref_crs)
teselas_vect <- vect(teselas_proj)

# ===============================
# 4. RASTERS & ALINEACIÓN DE CRS
# ===============================
r_files <- list.files(dir_aoh, pattern = "_AOH\\.tif$", full.names = TRUE)
if (length(r_files) == 0) stop("❌ No hay AOHs en la carpeta especificadas")

# Tomamos el CRS del primer ráster y transformamos las teselas a ese CRS (¡UNA SOLA VEZ!)
ref_crs <- crs(rast(r_files[1]))
teselas_proj <- st_transform(teselas, crs = ref_crs)
teselas_vect <- vect(teselas_proj)

# ===============================
# 5. LOOP PRINCIPAL
# ===============================
# Cambia '10' por 'nrow(teselas_vect)' para procesar todo el país
num_teselas_a_procesar <- nrow(teselas_vect)
resultados_lista <- vector("list", num_teselas_a_procesar)

cat("🚀 Iniciando cálculo STAR para", num_teselas_a_procesar, "teselas...\n")

for (i in seq_len(num_teselas_a_procesar)) {
  
  poly <- teselas_vect[i]
  poly_id <- teselas_proj$GRID_ID[i]
  ext_poly <- ext(poly)
  
  cat("🧩 Procesando Tesela", i, "/", num_teselas_a_procesar, "\n")
  resultados_poly <- list()
  
  for (j in seq_along(r_files)) {
    
    r_file <- r_files[j]
    spname_file <- basename(r_file)
    
    row_meta <- meta %>% filter(species_file == spname_file)
    if (nrow(row_meta) == 0) next
    
    peso <- row_meta$weight[1]
    if (is.na(peso) || peso == 0) next
    
    r <- rast(r_file)
    ext_r <- ext(r)
    
    # ⚡ FILTRO 1: Si las extensiones BBOX no intersecan, saltar rápido sin cargar ni recortar
    if (is.null(intersect(ext_poly, ext_r))) next
    
    # Convertir a binario si no lo está
    r_bin <- r == 1
    
    # ⚡ Recorte seguro
    r_crop <- tryCatch({ crop(r_bin, poly) }, error = function(e) NULL)
    if (is.null(r_crop)) next
    
    r_mask <- mask(r_crop, poly)
    
    # Conteo de celdas solapadas en la tesela
    ca_crop <- freq(r_mask)
    overlap_cells <- ifelse(any(ca_crop$value == 1, na.rm = TRUE),
                            ca_crop$count[ca_crop$value == 1], 0)
    
    if (overlap_cells == 0) next
    
    # Conteo global del hábitat de la especie
    ca_global <- freq(r_bin)
    global_cells <- ifelse(any(ca_global$value == 1, na.rm = TRUE),
                           ca_global$count[ca_global$value == 1], 0)
    
    if (global_cells == 0) next
    
    # Cálculos STAR
    START <- overlap_cells / global_cells
    STAR_species <- START * peso
    
    resultados_poly[[j]] <- data.frame(
      GRID_ID       = poly_id,
      species       = row_meta$Species[1],
      file          = spname_file,
      iucn_status   = row_meta$iucn_status[1],
      weight        = peso,
      global_cells  = global_cells,
      overlap_cells = overlap_cells,
      START         = START,
      STAR_species  = STAR_species
    )
  }
  
  resultados_lista[[i]] <- bind_rows(resultados_poly)
}

# ===============================
# 6. UNIR RESULTADOS Y LIMPIAR
# ===============================
resultados_final <- bind_rows(resultados_lista)

if (nrow(resultados_final) == 0) {
  warning("⚠️ No se encontraron solapamientos de especies en las teselas analizadas.")
  STAR_por_tesela <- data.frame(GRID_ID = teselas$GRID_ID, STAR_total = 0)
} else {
  STAR_por_tesela <- resultados_final %>%
    group_by(GRID_ID) %>%
    summarise(STAR_total = sum(STAR_species, na.rm = TRUE), .groups = "drop")
}

# ===============================
# 7. UNIR CON SHAPE (Manejando NAs)
# ===============================
teselas_STAR <- teselas %>%
  left_join(STAR_por_tesela, by = "GRID_ID") %>%
  mutate(STAR_total = ifelse(is.na(STAR_total), 0, STAR_total)) # Rellenar NAs con 0

# ===============================
# 8. GUARDAR
# ===============================
# Para guardar como Shapefile:
st_write(teselas_STAR, 
         "C:/Users/walter.garcia/Documents/Piloto STAR Colombia/star_t2.shp", 
         delete_dsn = TRUE)

write.csv(resultados_final, 
          "C:/Users/walter.garcia/Documents/Piloto STAR Colombia/STAR_teselas_detalle.csv", 
          row.names = FALSE)

cat("\n✅ STAR por teselas completado exitosamente.\n")



library(terra)
library(sf)
library(dplyr)
library(readr)
library(ggplot2)

# ===============================
# 1. RUTAS
# ===============================
dir_aoh  <- "C:/Users/walter.garcia/Documents/Piloto STAR Colombia/AOH"
csv_file <- "C:/Users/walter.garcia/Documents/STAR/Especies_amenazadas.csv"
aoi_file <- "C:/Users/walter.garcia/Documents/Piloto STAR Colombia/teselas_hexagonos_colombia_50km.shp"

# ===============================
# 2. METADATA & PESOS IUCN
# ===============================
meta <- read_delim(csv_file, delim = ";", show_col_types = FALSE)
meta$species_file <- basename(meta$species_file)

iucn_weights <- c("LC" = 0, "NT" = 100, "VU" = 200, "EN" = 300, "CR" = 400)
meta$weight <- iucn_weights[meta$iucn_status]

# ===============================
# 3. TESELAS
# ===============================
teselas <- st_read(aoi_file, quiet = TRUE)
teselas$GRID_ID <- 1:nrow(teselas)

# ===============================
# 4. RASTERS & CRS
# ===============================
r_files <- list.files(dir_aoh, pattern = "_AOH\\.tif$", full.names = TRUE)
if (length(r_files) == 0) stop("❌ No hay AOHs en la carpeta")

ref <- rast(r_files[1])

# Reproyectamos la geometría vectorial al CRS del ráster (¡Una sola vez!)
teselas_proj <- st_transform(teselas, crs(ref))
teselas_vect <- vect(teselas_proj)

# ===============================
# 5. PRECALCULAR ÁREA GLOBAL
# ===============================
cat("📦 Precalculando áreas globales de hábitat...\n")
global_stats <- list()

for (f in r_files) {
  r <- rast(f)
  r_bin <- ifel(r == 1, 1, NA)
  
  # Suma directa de celdas
  g_cells <- global(r_bin, "sum", na.rm = TRUE)[1, 1]
  global_stats[[basename(f)]] <- ifelse(is.na(g_cells), 0, g_cells)
}

# ===============================
# 6. LOOP PRINCIPAL (TODAS LAS TESELAS)
# ===============================
num_teselas <- nrow(teselas_vect)
resultados_lista <- vector("list", length = num_teselas)

cat("🚀 Iniciando cálculo STAR + Riqueza para", num_teselas, "teselas...\n")

for (i in seq_len(num_teselas)) {
  
  poly <- teselas_vect[i]
  poly_id <- teselas_proj$GRID_ID[i]
  ext_poly <- ext(poly)
  
  if (i %% 100 == 0 || i == 1) {
    cat("🧩 Progreso: Tesela", i, "/", num_teselas, "\n")
  }
  
  resultados_poly <- list()
  
  for (j in seq_along(r_files)) {
    r_file <- r_files[j]
    spname_file <- basename(r_file)
    
    row_meta <- meta %>% filter(species_file == spname_file)
    if (nrow(row_meta) == 0) next
    
    peso <- row_meta$weight[1]
    if (is.na(peso) || peso == 0) next
    
    global_cells <- global_stats[[spname_file]]
    if (is.null(global_cells) || global_cells == 0) next
    
    r <- rast(r_file)
    ext_r <- ext(r)
    
    # ⚡ FILTRO BBOX: Si no se solapan en extensión, saltar
    if (is.null(intersect(ext_poly, ext_r))) next
    
    r_bin <- ifel(r == 1, 1, NA)
    
    r_crop <- tryCatch({ crop(r_bin, poly) }, error = function(e) NULL)
    if (is.null(r_crop)) next
    
    r_mask <- mask(r_crop, poly)
    
    overlap_cells <- global(r_mask, "sum", na.rm = TRUE)[1, 1]
    overlap_cells <- ifelse(is.na(overlap_cells), 0, overlap_cells)
    
    if (overlap_cells == 0) next
    
    presence <- 1
    START <- overlap_cells / global_cells
    STAR_species <- START * peso
    
    resultados_poly[[j]] <- data.frame(
      GRID_ID       = poly_id,
      species       = row_meta$Species[1],
      file          = spname_file,
      iucn_status   = row_meta$iucn_status[1],
      weight        = peso,
      global_cells  = global_cells,
      overlap_cells = overlap_cells,
      presence      = presence,
      START         = START,
      STAR_species  = STAR_species
    )
  }
  
  resultados_lista[[i]] <- bind_rows(resultados_poly)
}

# ===============================
# 7. UNIR RESULTADOS Y MÉTRICAS
# ===============================
resultados_final <- bind_rows(resultados_lista)

if (nrow(resultados_final) == 0) {
  warning("⚠️ No se registraron solapamientos.")
  fusion <- data.frame(
    GRID_ID = teselas$GRID_ID,
    STAR_total = 0, richness = 0, richness_weighted = 0,
    STAR_norm = 0, richness_norm = 0, richness_weighted_norm = 0,
    STAR_RICH = 0
  )
} else {
  # Métricas
  STAR_por_tesela <- resultados_final %>%
    group_by(GRID_ID) %>%
    summarise(STAR_total = sum(STAR_species, na.rm = TRUE), .groups = "drop")
  
  riqueza_simple <- resultados_final %>%
    group_by(GRID_ID) %>%
    summarise(richness = sum(presence, na.rm = TRUE), .groups = "drop")
  
  riqueza_ponderada <- resultados_final %>%
    group_by(GRID_ID) %>%
    summarise(richness_weighted = sum(weight * presence, na.rm = TRUE), .groups = "drop")
  
  # Combinar
  fusion <- data.frame(GRID_ID = teselas$GRID_ID) %>%
    left_join(STAR_por_tesela, by = "GRID_ID") %>%
    left_join(riqueza_simple, by = "GRID_ID") %>%
    left_join(riqueza_ponderada, by = "GRID_ID") %>%
    mutate(
      STAR_total        = ifelse(is.na(STAR_total), 0, STAR_total),
      richness          = ifelse(is.na(richness), 0, richness),
      richness_weighted = ifelse(is.na(richness_weighted), 0, richness_weighted)
    )
  
  # Normalización
  max_star <- max(fusion$STAR_total, na.rm = TRUE)
  max_rich <- max(fusion$richness, na.rm = TRUE)
  max_rw   <- max(fusion$richness_weighted, na.rm = TRUE)
  
  fusion$STAR_norm              <- if (max_star > 0) fusion$STAR_total / max_star else 0
  fusion$richness_norm          <- if (max_rich > 0) fusion$richness / max_rich else 0
  fusion$richness_weighted_norm <- if (max_rw > 0) fusion$richness_weighted / max_rw else 0
  
  # Índice Híbrido (70% STAR + 30% Riqueza Ponderada)
  fusion$STAR_RICH <- (fusion$STAR_norm * 0.7) + (fusion$richness_weighted_norm * 0.3)
}

# ===============================
# 8. UNIR AL SHAPE & EXPORTAR
# ===============================
teselas_STAR <- teselas %>%
  left_join(fusion, by = "GRID_ID")

# Guardar GeoPackage (Recomendado para conservar nombres de columnas completos)
st_write(teselas_STAR, 
         "C:/Users/walter.garcia/Documents/Piloto STAR Colombia/star_ric.shp", 
         delete_dsn = TRUE, quiet = TRUE)

# Guardar también como Shapefile
st_write(teselas_STAR, 
         "C:/Users/walter.garcia/Documents/STAR/STAR_RICH_teselas.shp", 
         delete_dsn = TRUE, quiet = TRUE)

write.csv(resultados_final, 
          "C:/Users/walter.garcia/Documents/STAR/STAR_RICH_detalle.csv", 
          row.names = FALSE)

cat("\n✅ Cálculo STAR + Riqueza completado exitosamente.\n")

# =========================================
# MAPA GGPLOT2 (ADAPTADO)
# =========================================

# Ajustamos la variable a graficar (Richness Weighted o STAR_RICH)
# Nota: Si lees desde .shp, 'richness_weighted' se recorta a 'rchnss_w' o 'rch_wghtd'
var_mapa <- teselas_STAR$richness_weighted

# Crear la clasificación por cortes
teselas_STAR <- teselas_STAR %>%
  mutate(
    clase = cut(
      richness_weighted,
      breaks = c(-1, 0, 1600, 2200, 2600, 3600, max(richness_weighted, 4600, na.rm = TRUE)),
      include.lowest = TRUE,
      labels = c(
        "Sin especies (0)",
        "1 - 1600",
        "1600 - 2200",
        "2200 - 2600",
        "2600 - 3600",
        "> 3600"
      )
    )
  )

mapa <- ggplot(teselas_STAR) +
  geom_sf(
    aes(fill = clase),
    color = "black",
    linewidth = 0.05
  ) +
  scale_fill_manual(
    values = c(
      "Sin especies (0)" = "#f7f7f7",
      "1 - 1600"         = "#d80027",
      "1600 - 2200"     = "#f4a582",
      "2200 - 2600"     = "#cccccc",
      "2600 - 3600"     = "#888888",
      "> 3600"          = "#333333"
    ),
    drop = FALSE,
    name = NULL
  ) +
  labs(
    title = "Riqueza Ponderada por Tesela (Piloto Colombia)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 12, face = "bold"),
    legend.position = "left",
    legend.key.height = unit(0.8, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

print(mapa)

# Guardar imagen en alta resolución
ggsave(
  "C:/Users/walter.garcia/Documents/STAR/Mapa_riqueza_ponderada.png",
  plot = mapa,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)
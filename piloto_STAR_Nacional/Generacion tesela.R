library(sf)
library(dplyr)

# 1. Cargar el Shapefile de Colombia
colombia <- st_read("C:/Users/walter.garcia/Downloads/Colombia_Terrestre_Marino/colombia_magna_Maritimo.shp")


# 2. Convertir a proyección plana en metros (MAGNA-SIRGAS Origen Nacional)
colombia_flat <- st_transform(colombia, crs = 9377)

# 3. Crear la cuadrícula de HEXÁGONOS
# cellsize = 50000 es la distancia entre lados opuestos del hexágono (~50 km)
# square = FALSE le indica a R que genere geometría hexagonal
grid_hex <- st_make_grid(colombia_flat, cellsize = 50000, square = FALSE)

# 4. Convertir a objeto sf y asignar ID a cada hexágono
grid_hex_sf <- st_sf(geometry = grid_hex) %>%
  mutate(tile_id = row_number())

# 5. Filtrar solo los hexágonos que intersecan la masa terrestre de Colombia
teselas_hex_colombia <- st_filter(grid_hex_sf, colombia_flat)

# 6. Guardar en formato Shapefile de hexágonos
st_write(teselas_hex_colombia, "teselas_hexagonos_colombia_50km.shp", delete_layer = TRUE)




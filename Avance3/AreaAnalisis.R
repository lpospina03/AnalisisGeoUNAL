library(spBayes)
library(maps)
library(RANN)
library(gjam)
library(CARBayes)
library(CARBayesdata)
library(mgcv)
library(spdep)
library(DClusterm)
data(NY8)
library(dplyr)
library(sf)
library(spdep)
library(INLA)
library(Matrix)
library(ggplot2)
library(cowplot)
library(ggspatial)
library(tidyr)
library(sf)         # Leer y manejar datos espaciales
library(lme4)       # Modelos multinivel (glmer)
library(pscl)       # Pseudo-R² para modelos GLM
library(MuMIn)      # R² en modelos multinivel
library(ggplot2)    # Gráficas
library(ggspatial)  # Agregar flechas de norte y escalas en mapas
library(sjPlot)     # Visualizar efectos de modelos multinivel
library(dplyr)      # Manipulación de datos
library(pROC)       # Curva ROC (para modelos binarios)
library(broom)      # Convertir objetos de modelos a data frames limpias

# Leer y transformar el geojson
datos <- st_read("C:/Users/USUARIO/Documents/AnalisisGeoEspacial/Codigo/AreaAnalisisFinal/AreaAnalisisFinal.geojson")
datos <- st_transform(datos, crs = 9377)

# Coordenadas de centroides para la matriz de vecindad por distancia
coords <- st_coordinates(st_centroid(datos))
distancias <- spDists(coords, longlat = FALSE)
vecinos_250 <- dnearneigh(coords, 0, 250000)  # 250 km = 250000 m


# Verificar si hay polígonos sin vecinos
if (any(card(vecinos_250) == 0)) {
  cat("¡Alerta! Hay unidades sin vecinos dentro de 250 km.\n")
}

# Matriz de adyacencia binaria
matriz_binaria <- nb2mat(vecinos_250, style = "B")
grafo <- as(matriz_binaria, "Matrix")

# Lista de pesos para Moran
listw <- nb2listw(vecinos_250, style = "B")

# Crear identificador único por fila
datos$id <- 1:nrow(datos)


# Modelo base sin estructura espacial
m0 <- inla(
  Total_Eventos ~ 1 + Ptotal_Prom + Total_Estaciones + Elev_promedio,
  offset = log(Area_km2),
  data = as.data.frame(datos),
  family = "nbinomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)
datos$m0 <- m0$summary.fitted.values[1:nrow(datos), "mean"]

# Modelo ICAR (Besag)
m1 <- inla(
  Total_Eventos ~ 1 + Ptotal_Prom + Total_Estaciones + Elev_promedio +
    f(id, model = "besag", graph = grafo),
  offset = log(Area_km2),
  data = as.data.frame(datos),
  family = "nbinomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)
datos$m1_ICAR <- m1$summary.fitted.values[, "mean"]

# Modelo BYM
m2 <- inla(
  Total_Eventos ~ 1 + Ptotal_Prom + Total_Estaciones + Elev_promedio +
    f(id, model = "bym", graph = grafo),
  offset = log(Area_km2),
  data = as.data.frame(datos),
  family = "nbinomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)
datos$m2_BYM <- m2$summary.fitted.values[, "mean"]

# Modelo Leroux
Q <- Diagonal(nrow(grafo), rowSums(grafo)) - grafo
Cmatrix <- Diagonal(nrow(datos)) - Q
m3 <- inla(
  Total_Eventos ~ 1 + Ptotal_Prom + Total_Estaciones + Elev_promedio +
    f(id, model = "generic1", Cmatrix = Cmatrix),
  offset = log(Area_km2),
  data = as.data.frame(datos),
  family = "nbinomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)
datos$m3_Leroux <- m3$summary.fitted.values[, "mean"]


# Residuos de Pearson
datos$res_m0 <- (datos$Total_Eventos - datos$m0) / sqrt(datos$m0)
datos$res_m1 <- (datos$Total_Eventos - datos$m1_ICAR) / sqrt(datos$m1_ICAR)
datos$res_m2 <- (datos$Total_Eventos - datos$m2_BYM) / sqrt(datos$m2_BYM)
datos$res_m3 <- (datos$Total_Eventos - datos$m3_Leroux) / sqrt(datos$m3_Leroux)

# Tests de autocorrelación espacial de los residuos (Moran)
moran_m0 <- moran.mc(datos$res_m0, listw = listw, nsim = 999)
moran_m1 <- moran.mc(datos$res_m1, listw = listw, nsim = 999)
moran_m2 <- moran.mc(datos$res_m2, listw = listw, nsim = 999)
moran_m3 <- moran.mc(datos$res_m3, listw = listw, nsim = 999)

# Mostrar resultados
cat("Moran para modelo base (m0):      p =", moran_m0$p.value, "\n")
cat("Moran para modelo ICAR (m1):      p =", moran_m1$p.value, "\n")
cat("Moran para modelo BYM (m2):       p =", moran_m2$p.value, "\n")
cat("Moran para modelo Leroux (m3):    p =", moran_m3$p.value, "\n")


# Reorganizar datos para graficar
df_residuos <- datos %>%
  st_drop_geometry() %>%
  select(Cuenca, res_m0, res_m1, res_m2, res_m3) %>%
  pivot_longer(cols = starts_with("res_"), names_to = "modelo", values_to = "residual")

ggplot(df_residuos, aes(x = Cuenca, y = residual, fill = Cuenca)) +
  geom_boxplot(outlier.shape = 1) +
  facet_wrap(~modelo) +
  theme_bw() +
  labs(title = "Residuos de Pearson por modelo y cuenca")

cuenca_colores <- c("Cauca" = "red", "Magdalena" = "blue")

p_mapa <- ggplot() + 
  geom_sf(data = datos, aes(fill = res_m1, color = Cuenca), size = 0.6) +
  scale_fill_gradient2(low = "yellow", mid = "white", high = "purple", midpoint = 0, name = "Res. Pearson") +
  scale_color_manual(values = cuenca_colores) +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(1, "cm"), width = unit(1, "cm")) +
  theme_bw()

p_box <- ggplot(datos, aes(x = Cuenca, y = res_m1, fill = Cuenca)) +
  geom_boxplot(outlier.shape = 1) +
  scale_fill_manual(values = cuenca_colores) +
  theme_bw()

plot_grid(p_mapa, p_box, labels = c("A", "B"), ncol = 2)


#__________________________________________
### Modelo Jerarquico:

library(sf)
library(lme4)
library(dplyr)

# Leer el geojson
datos <- st_read("C:/Users/USUARIO/Documents/AnalisisGeoEspacial/Codigo/AreaAnalisisFinal/AreaAnalisisFinal.geojson")

# Asegurar que Cuenca sea factor
datos$Cuenca <- as.factor(datos$Cuenca)

# Estandarizar variables predictoras
datos <- datos %>%
  mutate(
    Ptotal_Prom_std = scale(Ptotal_Prom),
    Total_Estaciones_std = scale(Total_Estaciones),
    Elev_promedio_std = scale(Elev_promedio)
  )

#Intercepto aleatoreo por cuenca: 

m_intercept <- glmer.nb(
  Total_Eventos ~ Ptotal_Prom_std + Total_Estaciones_std + Elev_promedio_std +
    (1 | Cuenca),
  data = datos,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(m_intercept)

# Modelo con pendiente aleatoria solo para una variable ( Ptotal_Prom_std):
m_pendiente <- glmer.nb(
  Total_Eventos ~ Ptotal_Prom_std + Total_Estaciones_std + Elev_promedio_std +
    (0 + Ptotal_Prom_std | Cuenca),
  data = datos,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(m_pendiente)

# Modelo con intercepto y pendientes aleatorias por cuenca:
m_intercept_pendientes <- glmer.nb(
  Total_Eventos ~ Ptotal_Prom_std + Total_Estaciones_std + Elev_promedio_std +
    (1 + Ptotal_Prom_std + Total_Estaciones_std + Elev_promedio_std | Cuenca),
  data = datos,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(m_intercept_pendientes)


##________________
# implementación de modelos jerárquicos en INLA

library(INLA)
library(sf)
library(spdep)

# Leer los datos
datos <- st_read("C:/Users/USUARIO/Documents/AnalisisGeoEspacial/Codigo/AreaAnalisisFinal/AreaAnalisisFinal.geojson")

# Convertir Cuenca a factor y crear índices únicos para INLA
datos$Cuenca <- as.factor(datos$Cuenca)
datos$id_cuenca <- as.numeric(datos$Cuenca)

# Estandarizar variables
datos$Ptotal_std <- scale(datos$Ptotal_Prom)
datos$Estaciones_std <- scale(datos$Total_Estaciones)
datos$Elevacion_std <- scale(datos$Elev_promedio)

# Crear interacciones para pendientes aleatorias
datos$Cuenca.Ptotal <- datos$id_cuenca
datos$Cuenca.Elev <- datos$id_cuenca
datos$Cuenca.Estaciones <- datos$id_cuenca

## Intercepto aleatorio por cuenca: 


m_intercepto <- inla(
  Total_Eventos ~ 1 + Ptotal_std + Estaciones_std + Elevacion_std +
    f(id_cuenca, model = "iid"),
  offset = log(Area_km2),
  data = as.data.frame(datos),
  family = "nbinomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(m_intercepto)
datos$fitted_intercepto <- m_intercepto$summary.fitted.values$mean


# Modelo 2: intercepto + pendientes aleatorias por cuenca

m_intercepto_pendientes <- inla(
  Total_Eventos ~ 1 + Ptotal_std + Estaciones_std + Elevacion_std +
    f(id_cuenca, model = "iid") +  # intercepto aleatorio
    f(Cuenca.Ptotal, Ptotal_std, model = "iid", group = id_cuenca) +  # pendiente aleatoria
    f(Cuenca.Estaciones, Estaciones_std, model = "iid", group = id_cuenca) +
    f(Cuenca.Elev, Elevacion_std, model = "iid", group = id_cuenca),
  offset = log(Area_km2),
  data = as.data.frame(datos),
  family = "nbinomial",
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE)
)

summary(m_intercepto_pendientes)
datos$fitted_intercepto_pendientes <- m_intercepto_pendientes$summary.fitted.values$mean

#residuos de Pearson y autocorrelación espacial (si tienes matriz de vecindad):
# Crear matriz de vecinos (si tienes geometría)
datos_9377 <- st_transform(datos, crs = 9377)
coords <- st_centroid(st_geometry(datos_9377)) |> st_coordinates()

# Matriz de vecindad basada en distancia de 250 km
nb_250km <- dnearneigh(coords, 0, 250000)  # en metros
listw_250km <- nb2listw(nb_250km, style = "B")

# Residuos de Pearson para el modelo completo
res_pearson <- (datos$Total_Eventos - datos$fitted_intercepto_pendientes) / sqrt(datos$fitted_intercepto_pendientes)
datos$res_pearson <- res_pearson

# Moran's I sobre residuos
moran_test <- moran.mc(res_pearson, listw = listw_250km, nsim = 999, alternative = "greater")
print(moran_test)



















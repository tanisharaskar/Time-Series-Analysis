library(readr)
exo <- read_csv("E:/Projects/Exoplanet/Exoplanets.csv")
library(dplyr)

vars_needed <- c(
  "planet_name",
  "discovery_year",
  "discovery_method",
  "discovery_facility",
  "planet_radius_earth_radius",
  "insolation_flux_earth_flux",
  "equilibrium_temperature_k", 
  "eccentricity",
  "orbital_period_days",
  "stellar_age_gyr",
  "stellar_effective_temp_k",
  "stellar_metallicity_dex"
)

df_clean <- exo %>% select(all_of(vars_needed))
df_clean <- df_clean %>%
  mutate(
    planet_radius_earth_radius = as.numeric(planet_radius_earth_radius),
    insolation_flux_earth_flux = as.numeric(insolation_flux_earth_flux),
    equilibrium_temperature_k  = as.numeric(equilibrium_temperature_k),
    eccentricity               = as.numeric(eccentricity),
    orbital_period_days        = as.numeric(orbital_period_days),
    stellar_age_gyr            = as.numeric(stellar_age_gyr),
    stellar_effective_temp_k   = as.numeric(stellar_effective_temp_k),
    stellar_metallicity_dex    = as.numeric(stellar_metallicity_dex)
  ) %>%
  mutate(
    planet_radius_earth_radius = ifelse(is.na(planet_radius_earth_radius), median(planet_radius_earth_radius, na.rm = TRUE), planet_radius_earth_radius),
    insolation_flux_earth_flux = ifelse(is.na(insolation_flux_earth_flux), median(insolation_flux_earth_flux, na.rm = TRUE), insolation_flux_earth_flux),
    equilibrium_temperature_k  = ifelse(is.na(equilibrium_temperature_k),  median(equilibrium_temperature_k,  na.rm = TRUE), equilibrium_temperature_k),
    eccentricity               = ifelse(is.na(eccentricity),               median(eccentricity,               na.rm = TRUE), eccentricity),
    orbital_period_days        = ifelse(is.na(orbital_period_days),        median(orbital_period_days,        na.rm = TRUE), orbital_period_days),
    stellar_age_gyr            = ifelse(is.na(stellar_age_gyr),            median(stellar_age_gyr,            na.rm = TRUE), stellar_age_gyr),
    stellar_effective_temp_k   = ifelse(is.na(stellar_effective_temp_k),   median(stellar_effective_temp_k,   na.rm = TRUE), stellar_effective_temp_k),
    stellar_metallicity_dex    = ifelse(is.na(stellar_metallicity_dex),    median(stellar_metallicity_dex,    na.rm = TRUE), stellar_metallicity_dex)
  )
colSums(is.na(df_clean))
# S_size: triangular score based on planet radius (in Earth radii)
Rmin  <- 0.5
Rpeak <- 1.0
Rmax  <- 1.8

S_size <- function(Rp) {
  case_when(
    Rp <= Rmin | Rp >= Rmax ~ 0,
    Rp < Rpeak              ~ (Rp - Rmin) / (Rpeak - Rmin),
    Rp >= Rpeak             ~ (Rmax - Rp) / (Rmax - Rpeak)
  )
}

df_clean <- df_clean %>%
  mutate(planet_radius_earth_radius = as.numeric(planet_radius_earth_radius))
df_clean <- df_clean %>%
  mutate(s_size = S_size(planet_radius_earth_radius))
summary(df_clean$s_size)

# S_flux: triangular score based on insolation flux (relative to Earth)
Smin  <- 0.3
Speak <- 1.0
Smax  <- 2.0

S_flux <- function(S) {
  case_when(
    S <= Smin | S >= Smax ~ 0,
    S < Speak             ~ (S - Smin) / (Speak - Smin),
    S >= Speak            ~ (Smax - S) / (Smax - Speak)
  )
}

df_clean <- df_clean %>%
  mutate(insolation_flux_earth_flux = as.numeric(insolation_flux_earth_flux),
         s_flux = S_flux(insolation_flux_earth_flux))

summary(df_clean$s_flux)

# S_orbit: eccentricity score
e0   <- 0.1
emax <- 0.6

S_orbit <- function(e) {
  case_when(
    e <= e0   ~ 1,
    e <= emax ~ 1 - (e - e0) / (emax - e0),
    e > emax  ~ 0
  )
}

df_clean <- df_clean %>%
  mutate(s_orbit = S_orbit(eccentricity))

summary(df_clean$s_orbit)

# S_star: stellar environment score
tmin  <- 0.5
tpeak <- 4.5
tmax  <- 10.0

S_age <- function(t) {
  case_when(
    t <= tmin | t >= tmax ~ 0,
    t < tpeak             ~ (t - tmin) / (tpeak - tmin),
    t >= tpeak            ~ (tmax - t) / (tmax - tpeak)
  )
}

# Stellar effective temp score (peaks at Sun-like 5778K)
S_teff <- function(T) {
  case_when(
    T <= 3000 | T >= 7500 ~ 0,
    T < 5778              ~ (T - 3000) / (5778 - 3000),
    T >= 5778             ~ (7500 - T) / (7500 - 5778)
  )
}

# Stellar metallicity score (peaks at solar metallicity = 0)
S_met <- function(feh) {
  case_when(
    feh <= -1.0 | feh >= 0.5 ~ 0,
    feh < 0                  ~ (feh + 1.0) / (0 + 1.0),
    feh >= 0                 ~ (0.5 - feh) / (0.5 - 0)
  )
}

df_clean <- df_clean %>%
  mutate(
    s_age  = S_age(stellar_age_gyr),
    s_teff = S_teff(stellar_effective_temp_k),
    s_met  = S_met(stellar_metallicity_dex),
    s_star = (s_age * s_teff * s_met)^(1/3)
  )

summary(df_clean$s_star)

df_clean <- df_clean %>%
  mutate(H = (s_size * s_flux * s_orbit * s_star)^(1/4))

summary(df_clean$H)

threshold <- 0.3

df_clean %>%
  summarise(
    total_planets    = n(),
    high_H_count     = sum(H > threshold),
    high_H_fraction  = mean(H > threshold)
  )

library(ggplot2)

# yearly mean H and high-index counts
yearly_summary <- df_clean %>%
  mutate(discovery_year = as.numeric(discovery_year)) %>%
  group_by(discovery_year) %>%
  summarise(
    mean_H      = mean(H),
    high_H_count = sum(H > threshold),
    total        = n(),
    high_H_frac  = high_H_count / total
  )

# Plot 1: mean H over time
ggplot(yearly_summary, aes(x = discovery_year, y = mean_H)) +
  geom_line() +
  geom_point() +
  labs(title = "Mean Habitability Index by Discovery Year",
       x = "Discovery Year", y = "Mean H")

# Plot 2: fraction of high-H planets over time
ggplot(yearly_summary, aes(x = discovery_year, y = high_H_frac)) +
  geom_line() +
  geom_point() +
  labs(title = "Fraction of High-H Planets by Discovery Year",
       x = "Discovery Year", y = "Fraction H > 0.3")

method_summary <- df_clean %>%
  group_by(discovery_method) %>%
  summarise(
    mean_H       = mean(H),
    high_H_count = sum(H > threshold),
    total        = n(),
    high_H_frac  = high_H_count / total
  ) %>%
  arrange(desc(mean_H))

print(method_summary)

ggplot(method_summary, aes(x = reorder(discovery_method, mean_H), 
                           y = mean_H, fill = mean_H)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_flip() +
  labs(title = "Mean Habitability Index by Discovery Method",
       x = "Discovery Method", y = "Mean H")

facility_summary <- df_clean %>%
  group_by(discovery_facility) %>%
  summarise(
    mean_H = mean(H),
    total  = n()
  ) %>%
  filter(total > 10) %>%  # only facilities with enough planets
  arrange(desc(mean_H))

ggplot(facility_summary, aes(x = reorder(discovery_facility, mean_H), 
                             y = mean_H, fill = mean_H)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "darkred") +
  coord_flip() +
  labs(title = "Mean Habitability Index by Discovery Facility",
       x = "Facility", y = "Mean H")
view(ex)

location_data <- read.csv('NFL_Stadium_Locations.csv')

distance_matrix <- matrix(rep(NA, times=32^2), ncol=32)
rownames(distance_matrix) <- location_data$Team
colnames(distance_matrix) <- location_data$Team

deg2rad <- function(deg) {
   return(deg*pi/180)
}

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula
gcd <- function(long1, lat1, long2, lat2) {
   long1 = deg2rad(long1)
   lat1 = deg2rad(lat1)
   long2 = deg2rad(long2)
   lat2 = deg2rad(lat2)
   R <- 6371 # Earth mean radius [km]
   delta_long <- (long2 - long1)
   delta_lat <- (lat2 - lat1)
   a <- sin(delta_lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta_long/2)^2
   c <- 2 * asin(min(1, sqrt(a)))
   d = R * c
   return(d) # Distance in km
}

for (row_team in location_data$Team) {
   for (col_team in location_data$Team) {
      row_team_lat <- location_data$latitude[location_data$Team == row_team]
      row_team_long <- location_data$longitude[location_data$Team == row_team]
      col_team_lat <- location_data$latitude[location_data$Team == col_team]
      col_team_long <- location_data$longitude[location_data$Team == col_team]
      dist <- gcd(row_team_long, row_team_lat, col_team_long, col_team_lat)
      distance_matrix[row_team, col_team] <- dist
   }
}

write.csv(distance_matrix, file='distance_matrix.csv', row.names=T)

#===============================================================================
#                                    PCA                                       #
#===============================================================================
#-------------------------------------Iris--------------------------------------

# Charger les packages nécessaires
library(ggplot2)

# Charger le jeu de données iris
data(iris)

# Sélectionner les variables explicatives (les mesures de la fleur)
X <- iris[, 1:4]

# Standardisation des données
scaled_X <- scale(X)

# Calculer la PCA
pca_result <- prcomp(scaled_X, center = TRUE, scale. = TRUE)
pca_result

# Extraire les composantes principales pour chaque observation
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data

# Ajouter la colonne 'Species' pour la couleur
pca_data$Species <- iris$Species
pca_data

# Créer le graphique en utilisant ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  labs(title = "Analyse en Composantes Principales (PCA) sur iris",
       x = "Première Composante Principale (PC1)",
       y = "Deuxième Composante Principale (PC2)") +
  theme_minimal()

rm(iris,pca_data,pca_result,scaled_X,X)

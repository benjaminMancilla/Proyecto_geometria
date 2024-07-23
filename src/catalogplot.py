#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 03:17:50 2024

@author: benjamin
"""

import pandas as pd
import matplotlib.pyplot as plt

# Lee el archivo CSV
clusters_df = pd.read_csv('clusters.csv')

# Obtener el valor máximo de Y para invertir el eje
max_y = clusters_df['CentroidY'].max()

# Crear la figura y el eje con mayor tamaño y resolución
fig, ax = plt.subplots(figsize=(10, 8), dpi=200)

# Graficar los círculos para cada cluster
for _, row in clusters_df.iterrows():
    cluster_id = int(row['ClusterID'])
    centroid_x = row['CentroidX']
    centroid_y = row['CentroidY']
    radius = row['Radius']
    
    # Invertir la coordenada Y
    inverted_y = max_y - centroid_y + 50
    
    # Crear un círculo
    circle = plt.Circle((centroid_x, inverted_y), radius, color='dimgray', fill=False, linewidth=1)
    ax.add_patch(circle)
    
    # Añadir el ID del cluster
    ax.text(centroid_x, inverted_y, str(cluster_id), color='dimgray', fontsize=12, ha='center', va='center')

# Configurar los límites del gráfico
plt.xlim(clusters_df['CentroidX'].min() - 50, clusters_df['CentroidX'].max() + 50)
plt.ylim(0, max_y + 50)  # Ajustar el límite superior del eje Y
plt.gca().set_aspect('equal', adjustable='box')
xticks = ax.get_xticks()
print(xticks)
yticks = ax.get_yticks()
print(yticks)
xmax = 33/2
xmin = -33/2
ymax = 27/2
ymin = -27/2
xticks-=600
yticks-=450
xratio = 600/(33/2)
yratio = 450/(27/2)
print(xticks)
print(yticks)

def format_as_decimal(ticks,ratio):
    return [f'{tick/ratio:.{2}f}' for tick in ticks]

# Etiquetas de los ejes en unidades de minutos de arco
plt.xlabel('X Coordinate (arcminutes)')
plt.ylabel('Y Coordinate (arcminutes)')
ax.set_xticklabels(format_as_decimal(xticks, xratio))
ax.set_yticklabels(format_as_decimal(yticks, yratio))
# Título del gráfico
plt.title('Cluster Visualization')

# Mostrar el gráfico
plt.show()
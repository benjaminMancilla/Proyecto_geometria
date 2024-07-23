#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 14:57:10 2024

@author: benjamin
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
from scipy.interpolate import UnivariateSpline

# Leer los datos del archivo CSV
data = pd.read_csv('output.csv')

# Extraer las columnas necesarias
areas_normalizadas = data['Area/Promedio']
dp_values = data['dp(a)']

# Ordenar los datos
sorted_indices = np.argsort(areas_normalizadas)
sorted_areas = areas_normalizadas[sorted_indices]
sorted_dp_values = dp_values[sorted_indices]

# Crear una línea suave usando UnivariateSpline
spline = UnivariateSpline(sorted_areas, sorted_dp_values, s=1)  # Ajusta 's' para el nivel de suavizado
x_smooth = np.linspace(sorted_areas.min(), sorted_areas.max(), 500)
y_smooth = spline(x_smooth)

# Graficar puntos y línea suave
plt.figure(figsize=(10, 6))
plt.scatter(sorted_areas, sorted_dp_values, color='b', label='Datos', alpha=0.6)
plt.plot(x_smooth, y_smooth, color='r', label='Línea Suave')

# Agregar etiquetas y título
plt.xlabel('Área Normalizada')
plt.ylabel('Distribución dp(a)')
plt.title('Área Normalizada vs Distribución dp(a)')
plt.grid(True)
plt.legend()

# Mostrar la gráfica
plt.show()

# Calcular la densidad acumulativa
cum_density = np.arange(1, len(sorted_areas) + 1) / len(sorted_areas)

# Calcular la CDF de Poisson
mu = np.mean(sorted_areas)
poisson_cdf = poisson.cdf(sorted_areas, mu)

# Graficar densidad acumulativa y distribución de Poisson
plt.figure(figsize=(10, 6))
plt.plot(sorted_areas, cum_density, label='Densidad acumulativa', color='b')
plt.plot(sorted_areas, poisson_cdf, label='Distribución Poisson ajustada', linestyle='--', color='r')

# Añadir líneas verticales para los umbrales del 95% y 99%
plt.axvline(x=np.percentile(sorted_areas, 95), color='k', linestyle='-', label='Umbral 95%')
plt.axvline(x=np.percentile(sorted_areas, 99), color='k', linestyle='--', label='Umbral 99%')

# Agregar etiquetas y título
plt.xlabel('Área Normalizada')
plt.ylabel('Densidad acumulativa')
plt.title('Distribución de Densidad')
plt.legend()
plt.grid(True)

# Mostrar la gráfica
plt.show()



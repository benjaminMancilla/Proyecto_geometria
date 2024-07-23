# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 01:25:39 2024

@author: benja
"""

from astropy.io import fits
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
import csv

# Abre el archivo FITS
fits_file = '/home/benjamin/geo_proyect/src/CDS-P-ULTRAVISTA-color-KsJY~1.fits'
hdu_list = fits.open(fits_file)

# Iterar sobre cada canal
for i in range(hdu_list[0].data.shape[0]):
    image_data = hdu_list[0].data[i]

    # Estadísticas básicas de la imagen
    mean, median, std = sigma_clipped_stats(image_data, sigma=3.5)

    # Detectar fuentes en la imagen
    daofind = DAOStarFinder(fwhm=1.01, threshold=5.*std)
    sources = daofind(image_data - median)

    # Mostrar las posiciones de las fuentes detectadas
    print(f"Fuentes detectadas en el canal {i+1}:")
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # para una mejor visualización
        print(col, sources[col])

    # Mostrar la imagen y marcar las fuentes detectadas
    plt.figure(figsize=(10, 10))
    norm = simple_norm(image_data, 'sqrt', percent=99.5)
    plt.imshow(image_data, norm=norm, cmap='gray')
    plt.scatter(sources['xcentroid'], sources['ycentroid'], s=30, edgecolor='red', facecolor='none')
    plt.xlabel('X Pixel')
    plt.ylabel('Y Pixel')
    plt.title(f'Canal {i+1}')
    plt.show()
    
    csv_file = f'fuentes_detectadas_canal{i}.csv'

    # Abrir el archivo CSV en modo de escritura
    with open(csv_file, 'w', newline='') as file:
        writer = csv.writer(file)
    
        # Escribir encabezado opcional si se desea
        # writer.writerow(['xcentroid', 'ycentroid'])
    
        # Escribir las posiciones de las fuentes detectadas
        for source in sources:
            xcentroid = source['xcentroid']
            ycentroid = source['ycentroid']
            writer.writerow([xcentroid, ycentroid])
    
    print(f"Se han guardado las posiciones de las fuentes detectadas en '{csv_file}'.")

# Cierra el archivo FITS
hdu_list.close()
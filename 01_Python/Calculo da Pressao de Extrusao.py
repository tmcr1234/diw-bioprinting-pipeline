"""
Calculo da Pressao de Extrusao.py
=================================

╔══════════════════════════════════════════════════════════════════════════╗
║  ⚠ LEGACY (.xls) + ANTON PAAR MIGRATION PENDING                          ║
║  Built for legacy .xls layout (sheets named "Ramp*", stress at column    ║
║  index 2, viscosity at column index 4 after skiprows=3).                 ║
║  ALSO: file_path / save_path / samples below still point to an UNRELATED ║
║  legacy project. You MUST edit them before first use.                    ║
║  See Export/docs/01_rheology_data_formats.md                             ║
╚══════════════════════════════════════════════════════════════════════════╝

Quick Hagen-Poiseuille extrusion-pressure estimate from the peak viscosity
of a flow curve. Use as a sanity check, not for publication-grade numbers.
"""

import os
import pandas as pd
import xlrd
import numpy as np


# Define the dimensions
l_1_25_pol = 0.03175  # Comprimento da Agulha 1 1/4 Polegadas (m)
d_int_21_g = 0.000515  # Diametro Interno da agulha 21 Gauge (m)
l_1_pol = 0.0254  # Comprimento da Agulha 1 Polegadas (m)
d_int_22_g = 0.000413  # Diametro Interno da agulha 22 Gauge (m)
d_seringa = 0.0143 # Diametro Interno da Seringa BD 10ml (m)
l_seringa = 0.09 # Comprimento da Seringa BD 10 ml (m)
razao22G = l_1_pol/d_int_22_g
razao21G = l_1_25_pol/d_int_21_g
razaoSeringa = l_seringa/d_seringa

# Calculate the cross-sectional areas (A = pi * (d/2)^2)
# These areas are needed to calculate force (Force = Stress * Area)
area_21G = np.pi * (d_int_21_g / 2)**2
area_22G = np.pi * (d_int_22_g / 2)**2
area_seringa = np.pi * (d_seringa / 2)**2

# Sample names
samples = ['Quit+Gel', 'Gelatina', 'Quitosana', 'Gel de Cabelo 1', 'Gel de Cabelo 2']

# Paths
file_path = "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Artigo Carol/Dados do Artigo/Analises/Dados"
save_path = "/Users/thiagorodrigues/Library/CloudStorage/OneDrive-Pessoal/Artigo Carol/Dados do Artigo/Analises/Resultados"

# Define the output file path
output_file = "Pressão.txt"
output_path = os.path.join(save_path, output_file)

# Open the output file in write mode
with open(output_path, "w") as file:
    # Iterate over each sample
    for sample in samples:
        file_name = f"Visco - {sample}.xls"
        full_path = os.path.join(file_path, file_name)

        if os.path.exists(full_path):
            file.write(f"Processing {file_name}...\n")

            # Read Excel file
            xls = xlrd.open_workbook(full_path)

            # Iterate over each sheet
            for sheet_name in xls.sheet_names():
                if sheet_name.startswith("Ramp"):
                    df = pd.read_excel(xls, sheet_name, header=None, skiprows=3)  # skip first 3 rows
                    
                    # Check if the dataframe is not empty
                    if not df.empty:
                        # Assuming column index 1 contains the viscosity values and column index 2 contains stress
                        viscosidade_col = 4  # Update this index based on your data structure
                        stress_col = 2        # Update this index based on your data structure
                        
                        # Find the row with the maximum viscosity
                        max_viscosidade_row = df.iloc[:, viscosidade_col].idxmax()
                        
                        # Get the corresponding stress value from the row with max viscosity
                        stress = df.iloc[max_viscosidade_row, stress_col]
                        
                        # Calculate pressure
                        pressao_21G = 4 * razao21G * stress
                        pressao_22G = 4 * razao22G * stress
                        pressao_Seringa = 4 * razaoSeringa * stress

                        # Calculate the force (in Newtons) for each component
                        # Force (N) = Stress (Pa) * Area (m^2)
                        force_21G = stress * area_21G
                        force_22G = stress * area_22G
                        force_seringa = stress * area_seringa
                        
                        file.write(f"Sample: {sample}, Max Viscosidade Row: {max_viscosidade_row}, Stress: {stress} Pa, Pressão 21G: {pressao_21G} Pa, Pressão 22G: {pressao_22G} Pa, Pressão Seringa: {pressao_Seringa} Pa, Force 21G: {force_21G:.4f} N, Force 22G: {force_22G:.4f} N, Force Syringe: {force_seringa:.4f} N\n")
                    else:
                        file.write(f"No data found in sheet {sheet_name} for sample {sample}.\n")
        else:
            file.write(f"File not found: {file_name}\n")

print("Data saved to:", output_path)

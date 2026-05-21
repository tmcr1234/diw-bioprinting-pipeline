"""
Recovery_v1.py
==============

╔══════════════════════════════════════════════════════════════════════════╗
║  ⚠ LEGACY — superseded by extract_recovery_v2.py (Anton Paar 3iTT).      ║
║  This script computes recovery from a steady-shear ramp protocol         ║
║  (multiple shear rates, η-based). The new 3iTT-Osc-Rot-Osc protocol      ║
║  used with Anton Paar Rheocompass needs different math AND a different   ║
║  data layout (3 CSVs per sample). Use extract_recovery_v2 for any new    ║
║  data; keep this only to reproduce pre-migration numbers.                ║
║  See Export/docs/01_rheology_data_formats.md                             ║
╚══════════════════════════════════════════════════════════════════════════╝

Structural recovery % from the 3-step recovery test:
  Recovery (%) = viscosity(Ramp 3) / viscosity(Ramp 1) * 100
with first-order error propagation through the division.
"""

import os
import math
import pandas as pd

# =============================================================================
# CONFIGURAÇÕES
# =============================================================================
samples = ['C15', 'C15 4 NEG', 'Gel de Cabelo Bozano']
frequencies = [50, 100, 150, 200, 250, 300] 

# Caminhos dos arquivos
file_path = "./Reologia/Recovery"  
save_path = "./Analises/Python/Results"  

# Arquivos de saída
OUTPUT_CSV = os.path.join(save_path, "Dados_Recovery_Estruturado.csv")
OUTPUT_TXT = os.path.join(save_path, "Relatorio_Recovery.txt")

# =============================================================================
# EXTRAÇÃO DE DADOS E CÁLCULO
# =============================================================================
os.makedirs(save_path, exist_ok=True)

resultados = []
logs_de_erro = [] 

for sample in samples:
    for frequency in frequencies:
        file_name = f"{sample} - Recovery_{frequency}hz CP50 0,1mm.xls"
        full_path = os.path.join(file_path, file_name)

        if os.path.exists(full_path):
            try:
                xls = pd.ExcelFile(full_path, engine='xlrd')
                file_ramps = {} # Dicionário para organizar as abas deste arquivo
                
                # Lê todas as abas que começam com "Ramp"
                for sheet_name in xls.sheet_names:
                    if sheet_name.startswith("Ramp"):
                        df = pd.read_excel(xls, sheet_name=sheet_name, header=None, skiprows=3)
                        
                        if not df.empty and df.shape[1] > 4:
                            file_ramps[sheet_name] = {
                                'mean': df.iloc[:, 4].mean(),
                                'std': df.iloc[:, 4].std()
                            }
                        else:
                            logs_de_erro.append(f"{file_name} ({sheet_name}): Dados insuficientes.")
                
                # Garante a ordem correta das rampas (ex: Ramp 1, Ramp 2, Ramp 3)
                ramp_keys = sorted(list(file_ramps.keys()))
                
                # Precisamos de pelo menos a rampa inicial e a rampa final para o cálculo
                if len(ramp_keys) >= 3:
                    aba_inicial = ramp_keys[0] # (i) low shear, 25s
                    aba_meio = ramp_keys[1]    # (ii) high shear
                    aba_final = ramp_keys[2]   # (iii) recovery, 250s
                    
                    mean_1, std_1 = file_ramps[aba_inicial]['mean'], file_ramps[aba_inicial]['std']
                    mean_2, std_2 = file_ramps[aba_meio]['mean'], file_ramps[aba_meio]['std']
                    mean_3, std_3 = file_ramps[aba_final]['mean'], file_ramps[aba_final]['std']
                    
                    # Cálculo do Recovery (%)
                    if mean_1 > 0:
                        recovery_pct = (mean_3 / mean_1) * 100
                        # Propagação de erro da divisão
                        erro_pct = recovery_pct * math.sqrt((std_1/mean_1)**2 + (std_3/mean_3)**2)
                        
                        resultados.append({
                            "Amostra": sample,
                            "Shear_Rate_s-1": frequency,
                            "Visc_R1_Pa_s": round(mean_1, 4),
                            "Visc_R2_Pa_s": round(mean_2, 4),
                            "Visc_R3_Pa_s": round(mean_3, 4),
                            "Recovery_%": round(recovery_pct, 2),
                            "Erro_Recovery_%": round(erro_pct, 2)
                        })
                    else:
                        logs_de_erro.append(f"{file_name}: Viscosidade inicial (Ramp 1) igual a zero.")
                else:
                    logs_de_erro.append(f"{file_name}: Abas insuficientes para cálculo do Recovery (Encontrado: {len(ramp_keys)}).")
                    
            except Exception as e:
                logs_de_erro.append(f"{file_name}: Erro ao processar arquivo -> {e}")
        else:
            logs_de_erro.append(f"{file_name}: Arquivo não encontrado.")

# =============================================================================
# GERAÇÃO DO CSV
# =============================================================================
df_out = pd.DataFrame(resultados)
df_out.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")

# =============================================================================
# GERAÇÃO DO RELATÓRIO TXT
# =============================================================================
with open(OUTPUT_TXT, "w", encoding="utf-8") as f:
    f.write("=" * 115 + "\n")
    f.write("  RELATÓRIO DE EXTRAÇÃO DE DADOS -- TESTE DE RECUPERAÇÃO DE FORMA (NeG)\n")
    f.write("=" * 115 + "\n\n")

    f.write("MEMÓRIA DE EXTRAÇÃO E CÁLCULO\n")
    f.write("  - Arquivos Lidos: [Amostra] - Recovery_[Shear Rate]hz CP50 0,1mm.xls\n")
    f.write("  - Rampa Inicial (R1): Viscosidade média do primeiro intervalo (baixo cisalhamento).\n")
    f.write("  - Rampa Final (R3): Viscosidade média do último intervalo (recuperação).\n")
    f.write("  - Fórmula de Recovery: (Viscosidade_R3 / Viscosidade_R1) * 100\n")
    f.write("  - Propagação de Erro: Calculada através da raiz quadrada das incertezas relativas ao quadrado.\n\n")

    f.write("TABELA DE RESULTADOS\n")
    
    # Cabeçalho ajustado
    header = f"  {'Amostra':<23} {'Shear Rate':<12} {'Visc R1 (Pa.s)':<16} {'Visc R3 (Pa.s)':<16} {'Recovery (%)':<15} {'Erro (+/- %)':<15}\n"
    f.write(header)
    f.write("  " + "-" * 102 + "\n") 
    
    for r in resultados:
        f.write(
            f"  {r['Amostra']:<23} " 
            f"{r['Shear_Rate_s-1']:<12} " 
            f"{r['Visc_R1_Pa_s']:<16.4f} "
            f"{r['Visc_R3_Pa_s']:<16.4f} "
            f"{r['Recovery_%']:<15.2f} "
            f"{r['Erro_Recovery_%']:<15.2f}\n"
        )
    f.write("  " + "-" * 102 + "\n\n")

    if logs_de_erro:
        f.write("LOG DE AVISOS / ARQUIVOS AUSENTES\n")
        for erro in logs_de_erro:
            f.write(f"  - {erro}\n")
    else:
        f.write("LOG DE AVISOS\n")
        f.write("  - Todos os arquivos esperados foram encontrados, lidos e calculados com sucesso.\n")

print(f"✅ Processamento concluído!\nResultados salvos em:\n  -> {OUTPUT_CSV}\n  -> {OUTPUT_TXT}")
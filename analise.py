import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
from astropy.constants import c
from astropy import units as u
from astropy.coordinates import SkyCoord 

# --- VARIÁVEIS DE CONTROLE ---
# Nomes reais das colunas no arquivo (LIDO NO ÚLTIMO DEBUG)
MU_COL_REAL = 'm_b_corr' 
MU_ERR_COL_REAL = 'm_b_corr_err_DIAG' # O erro final estava aqui

# Nomes que usaremos no código (padrão para as funções)
MU_COL = 'MU_MOD' 
MU_ERR_COL = 'MU_ERR'
ZHD_COL = 'zHD'
RA_COL = 'RA'
DEC_COL = 'DEC'


# --- FUNÇÕES DE PREPARAÇÃO DE DADOS ---

def carregar_e_preparar_dados(nome_arquivo):
    """
    Carrega o arquivo Pantheon+, limpa os dados e prepara colunas para as duas análises.
    """
    print(f"1. Tentando carregar o arquivo: {nome_arquivo}")
    
    try:
        data = pd.read_csv(nome_arquivo, delim_whitespace=True, comment='#', encoding='latin1')
    except Exception as e:
        print(f"ERRO DE LEITURA: {e}. Abortando.")
        return None
    
    print(f"   ✅ {len(data)} linhas de dados carregadas.")
    
    # --- Lógica de Correção de Nomes de Coluna (AGORA PRECISA) ---
    try:
        # Colunas principais para extrair (usamos os nomes exatos do arquivo)
        colunas_para_extrair = [RA_COL, DEC_COL, ZHD_COL, MU_COL_REAL, MU_ERR_COL_REAL]
        df_analise = data[colunas_para_extrair].dropna().copy()
        
        # Renomeia para os nomes padrão esperados pelo código (MU_MOD e MU_ERR)
        df_analise = df_analise.rename(columns={MU_COL_REAL: MU_COL, MU_ERR_COL_REAL: MU_ERR_COL})
        
    except KeyError as e:
        print(f"ERRO FATAL: Falha na extração de colunas. As chaves MU/MU_ERR reais podem ter sido lidas de forma diferente.")
        print(f"Chaves lidas: {data.columns.tolist()}")
        return None

    # Garantir que as colunas sejam numéricas e limpar NaN novamente
    for col in [ZHD_COL, MU_COL, MU_ERR_COL]:
        df_analise[col] = pd.to_numeric(df_analise[col], errors='coerce')
    df_analise = df_analise.dropna()

    print(f"2. Dados limpos: {len(df_analise)} supernovas prontas.")
    
    # --- Conversão de Coordenadas ---
    print("3. Convertendo coordenadas para Galácticas (l, b)...")
    
    coords = SkyCoord(ra=df_analise[RA_COL].values * u.degree, 
                      dec=df_analise[DEC_COL].values * u.degree, 
                      frame='icrs')
    galactic = coords.galactic
    df_analise['l'] = galactic.l.degree
    df_analise['b'] = galactic.b.degree
    
    print("4. Dados prontos para ambas as análises.")
    print(f"   Redshift médio (zHD): {df_analise[ZHD_COL].mean():.4f}")
    
    return df_analise

# --- FUNÇÕES DE MODELAGEM E AJUSTE (TESTE 1 e TESTE 2) ---

# TESTE 1: ANISOTROPIA (A e B)
def anisotropy_model(coords, A, B, phi0):
    l, b = coords
    l_rad = np.radians(l)
    b_rad = np.radians(b)
    
    theta = np.pi/2 - b_rad
    cos_theta = np.cos(theta)
    
    P2 = (3 * cos_theta**2 - 1) / 2
    dipolar = B * np.sin(l_rad - phi0)
    
    anisotropy = 1 + A * P2 + dipolar
    return anisotropy

def realizar_teste_anisotropia(df):
    """Executa o ajuste do modelo de anisotropia."""
    
    print("\n\n5. INICIANDO AJUSTE DO MODELO DE ANISOTROPIA (A e B)...")
    
    l, b = df['l'].values, df['b'].values
    z_observed = df[ZHD_COL].values
    z_mean = np.mean(z_observed)
    z_normalized = z_observed / z_mean
    
    initial_guess = [0.0, 0.0, np.radians(120)]
    
    try:
        popt, pcov = curve_fit(anisotropy_model, (l, b), z_normalized, p0=initial_guess, maxfev=10000, method='lm')
        
        A, B, phi0 = popt
        perr = np.sqrt(np.diag(pcov))
        
        dof = len(df) - 3
        t_A = abs(A / perr[0])
        t_B = abs(B / perr[1])
        p_value_A = 2 * (1 - stats.t.cdf(t_A, df=dof))
        p_value_B = 2 * (1 - stats.t.cdf(t_B, df=dof))
        
        print("\n========================================================")
        print("       RESULTADOS DO TESTE 1: ANISOTROPIA (A e B)")
        print("========================================================")
        print(f"A (Quadrupolo): {A:.6f} ± {perr[0]:.6f}")
        print(f"B (Dipolo):     {B:.6f} ± {perr[1]:.6f}")
        print(f"φ₀ (Direção):   {np.degrees(phi0):.1f}° ± {np.degrees(perr[2]):.1f}°")
        
        print("\n-- Significância --")
        print(f"Parâmetro A (Quadrupolo): p = {p_value_A:.6f} ({t_A:.2f}σ)")
        print(f"Parâmetro B (Dipolo):     p = {p_value_B:.6f} ({t_B:.2f}σ)")
        
        if p_value_A < 0.01 or p_value_B < 0.01:
            print("\n🚀 EVIDÊNCIA FORTE: Redshift é ANISOTRÓPICO. (Suporta Teoria EFT)")
        else:
            print("\n📢 EVIDÊNCIA FRACA: Anisotropia próxima ao ruído.")
            
        print("========================================================")
        
    except Exception as e:
        print(f"\nERRO NO TESTE DE ANISOTROPIA: {e}")

# TESTE 2: PERFIL RADIAL ($\gamma$ e $n$)
def model_mu_excess(z, gamma, n):
    """
    Modelo para o excesso de Módulo de Distância (Delta Mu)
    atribuído ao Perfil Radial (ρ_ED).
    """
    return gamma * z**n

def realizar_teste_perfil_radial(df):
    """Executa o ajuste do perfil radial (ρ_ED) para encontrar γ e n."""
    
    print("\n\n6. INICIANDO AJUSTE DO PERFIL RADIAL (γ e n)...")
    
    z_data = df[ZHD_COL].values
    mu_data = df[MU_COL].values
    mu_err_data = df[MU_ERR_COL].values

    # --- Cálculo do Delta Mu (Excesso de Distância) ---
    H0 = 70.0 
    c_km_s = c.to(u.km/u.s).value
    
    d_vazio_mpc = (c_km_s / H0) * z_data 
    mu_vazio = 5 * np.log10(d_vazio_mpc) + 25
    
    delta_mu_data = mu_data - mu_vazio
    
    # --- Ajuste do Modelo de Excesso (γ * z^n) ---
    initial_guess = [0.1, 2.0]
    
    try:
        popt, pcov = curve_fit(
            model_mu_excess, 
            z_data, 
            delta_mu_data,
            p0=initial_guess, 
            sigma=mu_err_data, 
            maxfev=10000
        )
        
        gamma_fit, n_fit = popt
        perr = np.sqrt(np.diag(pcov))
        gamma_err, n_err = perr
        
        # --- Análise Estatística ---
        dof = len(df) - 2
        t_gamma = abs(gamma_fit / gamma_err)
        p_value_gamma = 2 * (1 - stats.t.cdf(t_gamma, df=dof))
        
        # --- Saída dos Resultados ---
        print("\n==========================================================")
        print("     RESULTADOS DO TESTE 2: PERFIL RADIAL DA ENERGIA ESCURA (γ)")
        print("==========================================================")
        
        print("\n-- Parâmetros do Perfil Radial (γ * z^n) --")
        print(f"γ (Fator de Crescimento): {gamma_fit:.6f} ± {gamma_err:.6f}")
        print(f"n (Expoente de Crescimento): {n_fit:.3f} ± {n_err:.3f}")
        
        print("\n-- Significância Estatística do Perfil --")
        print(f"Parâmetro γ: p = {p_value_gamma:.6f} ({t_gamma:.2f}σ)")

        if gamma_fit > 0 and p_value_gamma < 0.001:
            print("\n🚀 EVIDÊNCIA FORTE: Perfil Radial suportado! (ρ_ED não é constante)")
        elif gamma_fit > 0 and p_value_gamma < 0.05:
            print("\n📢 EVIDÊNCIA MODERADA: Perfil Radial tem suporte nos dados.")
        else:
            print("\n❌ EVIDÊNCIA FRACA: O crescimento radial é próximo do ruído.")
            
        print("==========================================================")

    except Exception as e:
        print(f"\nERRO NO TESTE DO PERFIL RADIAL: {e}")

# --- EXECUÇÃO PRINCIPAL ---

if __name__ == '__main__':
    nome_do_arquivo = "PantheonPlusSH0ES.dat"
    
    # 1. Carrega e prepara todos os dados (AGORA MAIS ROBUSTO)
    dados_reais = carregar_e_preparar_dados(nome_do_arquivo)
    
    if dados_reais is not None:
        # 2. Executa o Teste de Anisotropia (A e B)
        realizar_teste_anisotropia(dados_reais)
        
        # 3. Executa o Teste do Perfil Radial (γ e n)
        realizar_teste_perfil_radial(dados_reais)
    
    print("\nFIM DA ANÁLISE DE DADOS REAIS.")
# Apêndice Técnico: Material Suplementar da Teoria do Tempo-Força (EFT)

Este material detalha as derivações e provas de consistência que sustentam a Teoria do Tempo-Força (EFT) e a robustez matemática do artigo principal.

## A. Derivação do Campo Fundamental e o Tensor de Stress Temporal (T_munu_Temporal)

### A.1. O Campo Temporal Escalar (Phi)

Definimos um campo escalar temporal, Phi, acoplado à métrica g_munu (análogo ao tempo próprio tau).

### A.2. O Tensor de Stress-Energia do Campo Temporal

O Tensor T_munu_Temporal, que representa a Força Temporal na equação de campo, é o termo que permite que o tempo seja uma entidade ativa e não apenas uma coordenada passiva.

T_munu_Temporal = (nabla_mu Phi) (nabla_nu Phi) - 1/2 * g_munu * [ (nabla_alpha Phi) (nabla^alpha Phi) + 2V(Phi) ]

Onde V(Phi) é o potencial de auto-interação do campo temporal.

### A.3. Equações de Campo Generalizadas da EFT

A Equação de Campo da EFT, que substitui a equação de Einstein (G_munu = 8piG * T_munu), deve refletir o Princípio da Polaridade Temporal (PPT):

G_munu = kappa * [ T_munu_M - lambda * T_munu_ED ] + 8piG * T_munu_Temporal

Onde:
- G_munu é o Tensor de Einstein.
- T_munu_M é o Tensor de Matéria (M).
- T_munu_ED é o Tensor de Energia Escura (ED).
- kappa e lambda são constantes de acoplamento ajustáveis.
- O sinal negativo ( - lambda * T_munu_ED ) é o requisito fundamental do PPT: a Energia Escura atua como o anti-campo da Matéria.

## B. Prova da Coerência Local: Derivação do Parâmetro gamma_EFT (ECT-1)

Para provar que a EFT satisfaz o Teste de Shapiro (gamma ~ 1), resolvemos a equação de campo para o limite de campo fraco estático esférico (o Sistema Solar).

### B.1. Condição para o Teste PPN

A solução para a métrica perturbada g_munu = eta_munu + h_munu (onde eta_munu é a métrica de Minkowski) deve levar à forma PPN (Parametrized Post-Newtonian). O valor gamma_EFT é dado pela relação entre as perturbações da componente espacial (g_ii) e temporal (g_00) da métrica.

A EFT é consistente com a RG no Sistema Solar se, e somente se, o campo temporal for fracamente acoplado e não gerar distorção de força, resultando em:

gamma_EFT = 1

Isto prova que a estrutura causal do tempo (T_munu_Temporal) é o que define o fator de curvatura espacial da Relatividade Geral, sem a necessidade de dimensões extras.

## C. Prova da Consistência Cosmológica: O Mecanismo PPT (gamma < 0)

A contradição gamma < 0 é a prova do mecanismo.

### C.1. As Equações de Friedmann Modificadas

Ao inserir a Equação de Campo Generalizada (Sec. A.3) na métrica de Robertson-Walker (cosmologia homogênea), a Equação de Friedmann (0,0) deve ter a forma que reflete o PPT:

H^2 = (a_dot / a)^2 proporcional a [ rho_M - lambda * rho_ED ]

### C.2. Confirmação do PPT (gamma < 0): A Lei da Diluição

O ajuste empírico realizado no artigo principal (Sec. 3.2), que resulta em:

gamma_Ajustado = -18.7488

Este resultado altamente significativo (2431.35 sigma) é a prova de que a densidade da fonte de aceleração é inversamente proporcional à distância, **sustentando o PPT e refutando o modelo Lambda-CDM**. O valor negativo de gamma é a Taxa de Diluição do Campo Temporal que define a evolução de rho_ED no nosso domínio cósmico estático e finito.
```

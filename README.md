# silvdelta

Detecção e quantificação de mudanças estruturais em florestas a partir de
nuvens de pontos bitemporais (t0 = antes, t1 = depois da intervenção):
poda, desbaste e outras operações de manejo.

Funciona com nuvens TLS, MLS (SLAM) e ALS/UAV. O núcleo depende só de
`data.table`, `RANN` e `dbscan`; o `lidR` é opcional e usado apenas para ler
arquivos `.las/.laz`.

> Versão 0.2.0 — validada em talhões sintéticos (`simulate_stand()`).
> Ainda não validada com dados reais.

## Instalação

Direto do GitHub:

```r
install.packages("remotes")
remotes::install_github("nevesms/silvdelta")
```

Ou a partir do arquivo `.tar.gz`:

```r
install.packages(c("data.table", "RANN", "dbscan"))
install.packages("silvdelta_0.1.0.tar.gz", repos = NULL, type = "source")
# opcional, para ler LAS/LAZ:
install.packages("lidR")
```

## Uso rápido

```r
library(silvdelta)

# Totalmente automático: alvos de intensidade + operação identificada
res <- silvdelta("antes.las", "depois.las", registration = "targets")
res$operation        # desrama, desbaste (tipo), ambos, corte raso

# Árvore individual, com sementes (ex.: poda)
res <- silvdelta("antes.las", "depois.las",
                  seeds = "fustes.txt",   # X Y (Z) na ordem das árvores
                  level = "tree")
res
res$trees            # uma linha por árvore
plot(res)            # mapa
plot(res, "profile") # perfil vertical t0 x t1
export_results(res, "saida_silvdelta")

# Área, sem sementes (ex.: desbaste)
des <- silvdelta("antes.las", "depois.las", level = "area")
des$stand            # N/ha, G/ha, intensidade e tipo de desbaste
```

Sem dados à mão:

```r
sim <- simulate_stand(intervention = "both", seed = 1)
res <- silvdelta(sim, seeds = sim$seeds)
```

## Etapas do pipeline

| Etapa | Função | O que faz |
|---|---|---|
| 1. Épocas | `epochs()`, `read_epochs()` | Pareia as nuvens t0 e t1 |
| 2a. Registro por alvos | `detect_targets()`, `match_targets()`, `register_targets()`, `target_plot()` | Alvos retrorrefletivos achados pela intensidade; pareamento pelo arranjo; recorte, área e sistema local da parcela |
| 2b. Registro fino | `register_epochs()` | Alinha t1 em t0 usando só estruturas estáveis: ICP aparado na camada inferior dos fustes (posição, rotação, inclinação) e correção vertical pela diferença entre os terrenos |
| 3. Normalização | `normalize_epochs()` | Altura sobre um modelo de terreno **único** para as duas épocas |
| 4. Fustes | `detect_stems()`, `as_seeds()` | Localiza fustes pela continuidade vertical e ajusta círculo (RANSAC) para o DAP; com sementes ou automático |
| 5a. Árvore individual | `segment_epochs()`, `voxel_change()`, `tree_change()`, `classify_intervention()` | Voxels perdidos/ganhos com tolerância, atribuição por conectividade, altura de poda, classe por árvore |
| 5b. Área | `match_trees()`, `stand_change()`, `grid_change()`, `vertical_profile()` | Pareamento de fustes, resumo do povoamento, mudança em grade e perfil vertical |
| 6. Operação | `classify_operation()` | Nomeia a operação: desrama, desbaste (sistemático, por baixo, por cima, seletivo), ambos, corte raso |
| 7. Biomassa | `calibrate_biomass()` | Calibra volume perdido contra biomassa pesada em campo, com validação leave-one-out |

Cada etapa pode ser rodada separadamente; `silvdelta()` apenas encadeia todas.

## Alvos de intensidade (registro automático)

Instale um alvo retrorrefletivo em cada vértice da parcela e **um alvo extra
no meio de apenas uma das arestas**. Com `registration = "targets"`:

1. Pontos com intensidade acima de `mediana + 8 × MAD` (adaptado a cada
   levantamento e sensor) são agrupados; cada grupo compacto é um alvo.
   Reflexos isolados são descartados.
2. Os alvos de t0 e t1 são pareados pelas **distâncias entre eles**, que não
   dependem da posição nem da orientação de cada levantamento. As nuvens não
   precisam estar pré-alinhadas.
3. O alvo do meio da aresta quebra a simetria do retângulo, que sem ele
   encaixa igualmente bem girado 180°. Arranjos simétricos são sinalizados
   (`ambiguous`), com um aviso.
4. O ICP nas estruturas estáveis refina o alinhamento e corrige o Z pelo
   terreno.
5. Os vértices formam o polígono da parcela: recorte das duas épocas, área
   exata (em vez do envoltório da nuvem) e, com `plot_frame = TRUE`, um
   sistema de coordenadas local (origem num vértice, eixo X ao longo da
   aresta marcada).

Os pontos dos alvos e das hastes são removidos antes das análises. Se o seu
arquivo não tiver intensidade, informe os centros dos alvos manualmente em
`register_targets(targets_t0 = , targets_t1 = )`.

## Identificação da operação

`res$operation` traz a decisão e todas as métricas usadas:

| Resultado | Regra |
|---|---|
| `pruning` (desrama) | ≥ 30% das árvores em pé podadas (nível árvore); no nível área, `pruning_index` ≥ 0,3 (perda relativa na copa inferior menos a da copa superior) |
| `thinning` (desbaste) | ≥ 5% das árvores removidas |
| `thinning + pruning` | as duas condições |
| `clearcut` (corte raso) | ≥ 90% removidas |
| tipo `systematic` | ≥ 80% das removidas formam linhas de plantio completas |
| tipo `from below` / `from above` / `selective` | razão dg removidas / dg total < 0,9, > 1,1 ou entre os dois |

Todos os limiares são argumentos de `classify_operation()`.

## Sementes: com ou sem

- **Com sementes** (`seeds =` arquivo, data.frame ou matriz): cada semente é
  ligada ao fuste detectado mais próximo (raio `seed_radius`). A numeração das
  árvores segue as sementes, o que mantém o vínculo com o inventário de campo.
  Se nenhum fuste for encontrado perto da semente, a árvore permanece na
  análise com `found = FALSE`.
- **Sem sementes** (`seeds = NULL`): todos os fustes detectados viram árvores.

## Árvore individual x área

- `level = "tree"`: resultado por árvore (`pruned`, `removed`, `unchanged`),
  volume perdido, altura de poda, altura e base da copa em t0 e t1.
- `level = "area"`: fustes detectados de forma independente em t0 e t1 e
  pareados; resumo do povoamento (densidade, área basal, dg, intensidade do
  desbaste em N e em G, razão de desbaste dg removido / dg anterior).
- A grade (`$grid`), o perfil vertical (`$profile`) e o resumo do povoamento
  (`$stand`) são produzidos nos dois níveis.

## Decisões metodológicas

- **Registro em estruturas estáveis.** A copa muda com a intervenção; alinhar
  pela nuvem inteira enviesa o registro na direção da própria mudança.
- **Terreno comum.** Dois MDTs independentes criam diferença vertical falsa.
- **Tolerância na comparação de voxels.** Um voxel de t0 só é "perdido" se não
  houver voxel de t1 a menos de `tol` (padrão: maior entre 1 voxel e 2x o RMSE
  do registro). Grupos conectados pequenos viram "ruído".
- **Atribuição por conectividade.** Um grupo conectado de voxels perdidos (um
  galho) que toca um único fuste é atribuído inteiro àquela árvore, mesmo que a
  ponta invada o espaço da vizinha.
- **Altura de poda pelo padrão vertical.** A poda remove de baixo para cima;
  perdas espalhadas (oclusão) não formam essa faixa contínua e não são
  classificadas como poda.
- **Voxels "ganhos"** entre levantamentos próximos indicam principalmente
  oclusão em t0, e são reportados para avaliar essa assimetria.

## Limitações conhecidas (v0.1.0)

- Oclusão não é modelada (não usa a trajetória do sensor): material não visto
  em t1 conta como perdido.
- O eixo do fuste é tratado como vertical (árvores muito inclinadas perdem
  precisão nas métricas de copa).
- O modelo de terreno é simples (mínimos em grade); para terreno difícil,
  normalize antes e use `normalize = FALSE`.
- As nuvens precisam chegar aproximadamente alinhadas (erro < `max_dist`,
  padrão 0,5 m).
- `prune_height` é o alcance superior do material removido, que fica acima da
  altura de inserção do último galho podado quando os galhos são ascendentes.

## Próximos passos

1. Validar com dados reais de poda (inventário de campo) e de desbaste.
   Conferir a escala de intensidade do seu sensor nos alvos.
2. Máscara de oclusão a partir da trajetória/tempo GPS.
3. Eixo do fuste por fatias (inclinação e tortuosidade).
4. Separação madeira/folha no material removido.
5. Escrita de LAS com a classe de mudança por ponto; relatório automático.
6. Vinheta com estudo de caso; submissão ao CRAN/rOpenSci.

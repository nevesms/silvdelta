# silvdelta: fluxo passo a passo (equivalente a silvdelta(), com controle total)
# Troque o bloco "dados" pelos seus arquivos para rodar com dados reais.

library(silvdelta)
set.seed(1)

# ---- dados -----------------------------------------------------------------
sim   <- simulate_stand(nx = 4, ny = 3, intervention = "both", seed = 1)
ep    <- sim                # dados reais: ep <- read_epochs("antes.las", "depois.las")
seeds <- sim$seeds          # dados reais: seeds <- as_seeds("fustes.txt"), ou NULL

# ---- 1. registro em estruturas estaveis -------------------------------------
ep <- register_epochs(ep, dof = "rigid", layer = c(0.5, 3))
ep$registration$rmse        # erro residual (m) na camada estavel
ep$registration$dz          # correcao vertical aplicada (m)

# ---- 2. normalizacao com terreno comum --------------------------------------
ep <- normalize_epochs(ep, res = 0.5)

# ---- 3. fustes (com ou sem sementes) ----------------------------------------
stems <- detect_stems(ep$t0, seeds = seeds, tol = 0.02)
stems

# ---- 4a. arvore individual ---------------------------------------------------
ep    <- segment_epochs(ep, stems)
vc    <- voxel_change(ep, res = 0.05, stems = stems)
vc
trees <- classify_intervention(tree_change(ep, stems, vc))
trees[, .(tree_id, dbh, present_t1, lost_volume, prune_height, change)]
change_profile(vc, bin = 0.5)            # volume perdido por estrato

# ---- 4b. area ----------------------------------------------------------------
stems_t1 <- detect_stems(ep$t1)
pares    <- match_trees(stems, stems_t1)
stand_change(pares, area_ha = plot_area(ep$t0))
grid_change(ep, res = 1)
vertical_profile(ep, bin = 0.5)

# ---- 5. biomassa removida (com pesagem de campo) -----------------------------
# campo <- data.table::fread("inventario.csv")        # tree_id, biomassa_kg
# cal <- calibrate_biomass(trees$lost_volume, campo$biomassa_kg)
# cal                                                  # RMSE leave-one-out

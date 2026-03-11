# Quick AUC explainer figure
library(tidyverse)
library(pROC)
library(patchwork)
library(here)

# Load data
data_calanda = readRDS(here("Calanda_JSDM", "output", "data_calanda_jsdm_2026-03-06.rds"))
Y = data_calanda$Y

# Load out-of-fold predictions
fold_files = sort(list.files(
  here("Calanda_JSDM", "output", "results", "cv"),
  pattern = "^fold_.*a1_l0.01_le0.001\\.rds$", full.names = TRUE
))
P_oof = matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
for (f in fold_files) {
  fd = readRDS(f)
  P_oof[fd$test_idx, ] = fd$preds_test
}

# Pick a species near median AUC for illustration
species_aucs = sapply(1:ncol(Y), function(s) {
  tryCatch(as.numeric(auc(roc(Y[, s], P_oof[, s], quiet = TRUE))), error = function(e) NA)
})
median_auc = median(species_aucs, na.rm = TRUE)
sp_idx = which.min(abs(species_aucs - median_auc))
sp_name = colnames(Y)[sp_idx]
sp_auc = species_aucs[sp_idx]

y_true = Y[, sp_idx]
y_pred = P_oof[, sp_idx]

cat("Example species:", sp_name, "\n")
cat("AUC:", round(sp_auc, 3), "\n")
cat("Presences:", sum(y_true == 1), "/ Absences:", sum(y_true == 0), "\n")

# --- Panel A: Predicted probabilities by true class ---
df_pred = tibble(
  observed = factor(y_true, levels = c(0, 1), labels = c("Absent", "Present")),
  predicted = y_pred
)

p_a = ggplot(df_pred, aes(x = observed, y = predicted, fill = observed)) +
  geom_jitter(aes(color = observed), width = 0.25, alpha = 0.4, size = 1.2) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  scale_fill_manual(values = c("Absent" = "#81caf3", "Present" = "#d00000")) +
  scale_color_manual(values = c("Absent" = "#81caf3", "Present" = "#d00000")) +
  labs(
    title = "A) Predicted probabilities by true class",
    subtitle = "Good model: Present sites get higher predicted probabilities",
    x = "Observed", y = "Predicted probability"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# --- Panel B: Threshold sweep ---
thresholds = seq(0, 1, by = 0.01)
tpr_fpr = tibble(threshold = thresholds) %>%
  rowwise() %>%
  mutate(
    predicted_pos = sum(y_pred >= threshold),
    TP = sum(y_pred >= threshold & y_true == 1),
    FP = sum(y_pred >= threshold & y_true == 0),
    FN = sum(y_pred < threshold & y_true == 1),
    TN = sum(y_pred < threshold & y_true == 0),
    TPR = TP / (TP + FN),   # sensitivity / recall
    FPR = FP / (FP + TN)    # 1 - specificity
  ) %>%
  ungroup()

# Show 3 example thresholds
example_thresh = c(0.2, 0.5, 0.8)
df_examples = tpr_fpr %>% filter(threshold %in% example_thresh)

p_b = ggplot(tpr_fpr, aes(x = threshold)) +
  geom_line(aes(y = TPR, color = "TPR (Sensitivity)"), linewidth = 0.8) +
  geom_line(aes(y = FPR, color = "FPR (1 - Specificity)"), linewidth = 0.8) +
  geom_vline(data = tibble(threshold = example_thresh),
             aes(xintercept = threshold), linetype = "dotted", color = "grey40") +
  scale_color_manual(values = c("TPR (Sensitivity)" = "#d00000", "FPR (1 - Specificity)" = "#81caf3")) +
  labs(
    title = "B) TPR and FPR across thresholds",
    subtitle = "At each threshold: classify sites with P >= threshold as 'present'",
    x = "Classification threshold", y = "Rate", color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- Panel C: ROC curve ---
roc_obj = roc(y_true, y_pred, quiet = TRUE)
df_roc = tibble(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_c = ggplot(df_roc, aes(x = FPR, y = TPR)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = FPR, ymax = TPR), fill = "#d00000", alpha = 0.15) +
  geom_line(color = "#d00000", linewidth = 1) +
  geom_point(data = df_examples, aes(x = FPR, y = TPR), size = 3, color = "black") +
  geom_text(data = df_examples, aes(x = FPR, y = TPR,
            label = paste0("t=", threshold)),
            vjust = -1, size = 3.2) +
  annotate("text", x = 0.6, y = 0.3,
           label = paste0("AUC = ", round(sp_auc, 3)),
           size = 5, fontface = "bold", color = "#d00000") +
  annotate("text", x = 0.6, y = 0.2,
           label = "= shaded area\n= P(rank present > rank absent)",
           size = 3.2, color = "grey30") +
  labs(
    title = "C) ROC curve",
    subtitle = "Each point = one threshold; AUC = area under this curve",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw(base_size = 11) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

# --- Panel D: Intuition ---
# Random pairs: present vs absent
set.seed(42)
present_preds = y_pred[y_true == 1]
absent_preds  = y_pred[y_true == 0]
n_pairs = min(200, length(present_preds) * length(absent_preds))
pair_present = sample(present_preds, n_pairs, replace = TRUE)
pair_absent  = sample(absent_preds, n_pairs, replace = TRUE)
concordant = sum(pair_present > pair_absent)
tied = sum(pair_present == pair_absent)

df_pairs = tibble(
  pair = 1:n_pairs,
  present_pred = pair_present,
  absent_pred = pair_absent,
  result = case_when(
    present_pred > absent_pred ~ "Concordant",
    present_pred < absent_pred ~ "Discordant",
    TRUE ~ "Tied"
  )
)

p_d = ggplot(df_pairs, aes(x = absent_pred, y = present_pred, color = result)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("Concordant" = "#00bd89", "Discordant" = "#d00000", "Tied" = "grey50")) +
  annotate("text", x = 0.8, y = 0.2,
           label = paste0("Concordant: ", sum(df_pairs$result == "Concordant"), "/", n_pairs,
                          "\n= ", round(sum(df_pairs$result == "Concordant") / n_pairs, 3)),
           size = 3.5, color = "#00bd89", fontface = "bold") +
  labs(
    title = "D) AUC as pairwise comparison",
    subtitle = "200 random present-absent pairs: does present get higher P?",
    x = "Predicted P (absent site)", y = "Predicted P (present site)",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

# --- Combine ---
p_all = (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = paste0("How AUC works — example: ", sp_name, " (AUC = ", round(sp_auc, 3), ")"),
    subtitle = paste0("Out-of-fold predictions from 10-fold CV | ",
                      sum(y_true == 1), " presences, ", sum(y_true == 0), " absences"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    )
  )

pdf(here("Calanda_JSDM", "plot", "exp4_auc_explainer.pdf"), width = 14, height = 12)
print(p_all)
dev.off()
cat("Saved exp4_auc_explainer.pdf\n")

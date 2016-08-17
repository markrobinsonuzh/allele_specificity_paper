# Plot ROC curves for pure and noisy simulated data, comparing amrfinder, allelicmeth, and the ASM score.

# Pure data

load("../data/allelicmeth_amr_asm_roc_plot_pure.rda")

# plot ROC curves

pdf(file = "../figures/roc_amrfinder_allelicmeth_asm_score.pdf", w=10, h=8)

pred_amr <- prediction(roc_plot_list_2[[1]]$amr_score, roc_plot_list_2[[1]]$label)
pred_allelicmeth <- prediction(roc_plot_list_2[[1]]$allelicmeth_score, roc_plot_list_2[[1]]$label)
pred_asm_score <- prediction(roc_plot_list_2[[1]]$asm_score, roc_plot_list_2[[1]]$label)

perf_amr <- performance(pred_amr,  measure="sens", x.measure="fpr")
perf_allelicmeth <- performance(pred_allelicmeth,  measure="sens", x.measure="fpr")
perf_asm_score <- performance(pred_asm_score,  measure="sens", x.measure="fpr")

plot(perf_amr, col='darkgreen', lwd=2)
plot(perf_allelicmeth, add = TRUE, col='darkred', lwd=2)
plot(perf_asm_score, add = TRUE, col='darkblue', lwd=2)

for (i in 2:3) {
  pred_amr <- prediction(roc_plot_list_2[[i]]$amr_score, roc_plot_list_2[[i]]$label)
  pred_allelicmeth <- prediction(roc_plot_list_2[[i]]$allelicmeth_score, roc_plot_list_2[[i]]$label)
  pred_asm_score <- prediction(roc_plot_list_2[[i]]$asm_score, roc_plot_list_2[[i]]$label)
  
  perf_amr <- performance(pred_amr,  measure="sens", x.measure="fpr")
  perf_allelicmeth <- performance(pred_allelicmeth,  measure="sens", x.measure="fpr")
  perf_asm_score <- performance(pred_asm_score,  measure="sens", x.measure="fpr")
  
  plot(perf_amr, add = TRUE, col='darkgreen', lwd=2)
  plot(perf_allelicmeth, add = TRUE, col='darkred', lwd=2)
  plot(perf_asm_score, add = TRUE, col='darkblue', lwd=2)
}

abline(0,1,col='black')
legend(0.6,0.4,c('amrfinder', 'allelicmeth', 'ASM score'), lty=c(1,1), col=c('darkgreen', 'darkred', 'darkblue'))

dev.off()



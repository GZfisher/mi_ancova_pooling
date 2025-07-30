rm(list = ls())

################## compare the ANCOVA + pooling results ##################
load("ancova.RData")
r_res <- res_com100 %>% 
  mutate(AVISIT = paste0("WEEK",visit),
         AVISITN = visit,
         TRT01PN = 1,
         NImpute = n,
         Estimate = estimate,
         StdErr = SE,
         LCLMean = lower.CL,
         UCLMean = upper.CL,
         DF = df,
         tValue = t.ratio,
         Probt = p.value) %>% 
  select(AVISIT, AVISITN, TRT01PN, NImpute, Estimate, StdErr, LCLMean, UCLMean, DF, tValue, Probt, source)
sas_res <- haven::read_sas("SAS result/diff.sas7bdat") %>% 
  mutate(AVISITN = as.numeric(gsub("WEEK", "", AVISIT)),
         source = "SAS") %>% 
  arrange(AVISITN)
  
diffdf::diffdf(r_res, sas_res)

# haven::write_xpt(r_res, "diff_r.xpt", version = 5)

################## compare the MI results ##################
library(dplyr)
load("imputation.RData")
compare_sas <- haven::read_sas("imp_long.sas7bdat") %>% 
  select(SUBJID, TRT01PN, TRT01P, AVISIT, impno, CHG) %>% 
  mutate(AVISITN = as.numeric(gsub("WEEK", "", AVISIT))) %>% 
  left_join(data_complete %>% select(SUBJID, TRT01PN, TRT01P, AVISIT, AVISITN, CHG), 
            by = c("SUBJID", "TRT01P", "AVISIT", "AVISITN")) %>% 
  mutate(diff = CHG.x - CHG.y) %>% 
  filter(diff != 0) %>% 
  group_by(impno, AVISIT, AVISITN) %>% 
  summarise(mean_abs = mean(abs(diff)), mean_sqr = mean(diff**2), mean_r = mean(abs(diff/CHG.y))) %>% 
  group_by(AVISIT, AVISITN) %>% 
  summarise(mae = mean(mean_abs), mse = mean(mean_sqr), mrae = mean(mean_r)) %>% 
  arrange(AVISITN) %>% mutate(source = "SAS")


compare_r <- complete_long %>% filter(impno!=0) %>% 
  select(SUBJID, TRT01PN, TRT01P, AVISIT, AVISITN, impno, CHG) %>% 
  mutate(AVISITN = as.numeric(AVISITN)) %>% 
  left_join(data_complete %>% select(SUBJID, TRT01PN, TRT01P, AVISIT, AVISITN, CHG), 
        by = c("SUBJID", "TRT01P", "AVISIT", "AVISITN")) %>% 
  mutate(diff = CHG.x - CHG.y) %>% 
  filter(diff != 0) %>% 
  group_by(impno, AVISIT, AVISITN) %>% 
  summarise(mean_abs = mean(abs(diff)), mean_sqr = mean(diff**2), mean_r = mean(abs(diff/CHG.y))) %>% 
  group_by(AVISIT, AVISITN) %>% 
  summarise(mae = mean(mean_abs), mse = mean(mean_sqr), mrae = mean(mean_r)) %>% 
  arrange(AVISITN) %>% mutate(source = "R")

compare_all <- rbind(compare_sas, compare_r)

library(ggplot2)
p1 <- ggplot(compare_all, aes(x = AVISITN, y = mae, color = source, group = source)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "AVISIT",
       y = "MAE",
       color = "Source") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(compare_all, aes(x = AVISITN, y = mse, color = source, group = source)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "AVISIT",
       y = "MSE",
       color = "Source") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plot/MAE.png", p1, dpi = 300, width = 8, height = 6, units = "in")
ggsave("plot/MSE.png", p2, dpi = 300, width = 8, height = 6, units = "in")


################## compare the MI+ANCOVA results ##################

load("diff_r.RData")
diff_r1 <- diff_r %>% 
  mutate(AVISIT = paste0("WEEK",visit),
         AVISITN = visit,
         TRT01PN = 1,
         NImpute = n,
         Estimate = estimate,
         StdErr = SE,
         LCLMean = lower.CL,
         UCLMean = upper.CL,
         DF = df,
         tValue = t.ratio,
         Probt = p.value) %>% 
  select(AVISIT, AVISITN, TRT01PN, NImpute, Estimate, StdErr, LCLMean, UCLMean, DF, tValue, Probt, source)
diff_sas <- haven::read_sas("SAS result/diff_sas.sas7bdat") %>% 
  mutate(AVISITN = as.numeric(gsub("WEEK", "", AVISIT)),
         source = "SAS") %>% 
  arrange(AVISITN) %>% 
  select(AVISIT, AVISITN, TRT01PN, NImpute, Estimate, StdErr, LCLMean, UCLMean, DF, tValue, Probt, source)
diff_all <- rbind(diff_r1, diff_sas) %>% 
  mutate(AVISIT = factor(AVISIT, 
                         levels = unique(AVISIT)[order(as.numeric(gsub("WEEK", "", unique(AVISIT))))]))

# haven::write_xpt(diff_r1, "diff_r1.xpt", version = 5)

dotCOLS = c("#a6d8f0","#f9b282")
barCOLS = c("#008fd5","#de6b35")

# annotation <- data.frame(
#   y=rep(0.4,28), x=factor(seq(2,56,2)),levels = seq(2,56,2), 
#   label=paste0(mcmc_missingn[2:29,1],"  ",mono_missingn[2:29,1],"  ",all_missingn[2:29,1])
# )

p <- ggplot(diff_all, aes(x=AVISIT, y=Estimate)) + 
  #specify position here
  geom_linerange(linewidth=1,position=position_dodge(width = 0.5), aes(ymin=`LCLMean`, ymax=`UCLMean`,col=source)) +
  #geom_hline(yintercept=-0.8, lty=2) +
  #specify position here too
  geom_point(size=2, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5), aes(fill=source)) +
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  scale_x_discrete(name="AVISIT") +
  # geom_text(data=annotation, size=2, aes(x=x,y=y,label=label)) +
  # scale_y_continuous(name="Odds ratio", limits = c(0.5, 5)) +
  coord_flip() +
  theme_minimal()
p 
ggsave("plot/ANCOVA_forest.png", p, dpi = 300, width = 8, height = 4, units = "in")

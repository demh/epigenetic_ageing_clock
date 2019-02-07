###########################################################################################
#########                                                                         #########
#########                     Daniel Elias Martin Herranz                         #########
#########                             22/05/2018                                  #########
#########                              EMBL-EBI                                   #########
#########                           Thornton group                                #########
#########                                                                         #########
###########################################################################################

###########################################################################################
#####              Biological insights into the epigenetic ageing clock           #########
###########################################################################################
##### Comparison for the different reference-based cell-type deconvolution strategies. ####
##### Plot the results for the different predictions and decide on the best strategy.  ####
###########################################################################################
##### USAGE: manual                                                                    ####
###########################################################################################

###########################################################
##################### Dependencies ########################
###########################################################

library(Metrics);
library(ggplot2);
library(RColorBrewer);
library(gridExtra);


############################################################
################## Functions ###############################
############################################################

#### Function: assess the accuracy of the different methods to estimate cell-type composition. For each cell type in each method (as specified in Name column) it returns
#              the RMSE, MAE and R2.

# pred_df: dataframe that contains the different predicted values for cell-type composition (e.g. all_predictions). Columns: Name, Sample, Gran, CD4T, CD8T, B, Mono, NK. Data: percentages
# real_df: dataframe that contains the real values for cell-type composition (e.g. real_cc_gs). Columns: Sample, Gran, CD4T, CD8T, B, Mono, NK. Data: percentages

calculate_pred_error <- function(pred_df, real_df){
  
  for(m in unique(pred_df$Name)){ # For each method
    
    sub_method <- pred_df[pred_df$Name==m,];
    real_df_ordered <- real_df[match(sub_method$Sample, real_df$Sample),];
    rmses <- c();
    maes <- c();
    R2s <- c();
    
    for(i in 2:(ncol(real_df_ordered))){
      
      rmses <- c(rmses, rmse(actual=real_df_ordered[,i], predicted=sub_method[,i+1])); 
      maes <- c(maes, mae(actual=real_df_ordered[,i], predicted=sub_method[,i+1])); 
      R2s <- c(R2s, cor(real_df_ordered[,i], sub_method[,i+1])^2);
      
    }
    
    if(m==unique(pred_df$Name)[1]){ # First time, create final dataframe
      
      results_final <- data.frame(Name=rep(m, ncol(real_df_ordered)-1),
                                  Cell=colnames(real_df_ordered)[-1],
                                  RMSE=rmses, MAE=maes, R2=R2s);  
      
    }else{  # Append in final dataframe
      
      new_results <- data.frame(Name=rep(m, ncol(real_df_ordered)-1),
                                Cell=colnames(real_df_ordered)[-1],
                                RMSE=rmses, MAE=maes, R2=R2s); 
      results_final <- rbind(results_final, new_results);
    }
  }
  
  return(results_final);
}


############################################################
################## Running the pipeline ####################
############################################################

#### 1. Read the predictions and the actual values for the different samples in the goldstandard. 

all_predictions <- read.csv('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/compare_cc_strategies_predictions.csv');
real_cc_gs <- read.csv('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/compare_cc_strategies_real_values.csv');


#### 2. Estimate the performance of the different methods.

results_methods <- calculate_pred_error(all_predictions, real_cc_gs);


#### 3. Plot the results. 

jplot_RMSE <- ggplot(data=results_methods, aes(x=Name, y=RMSE, group=1)) + 
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean, geom="crossbar", size=0.5, color="grey") +
  geom_jitter(aes(col=Cell), position = position_jitter(height = 0, width = 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = 90, hjust = 1),
        axis.title=element_text(size=14,face="bold")) +
  scale_colour_brewer(type = "qual", palette = 'Dark2') + 
  xlab('Cell-type deconvolution strategy') + 
  geom_hline(yintercept=min(aggregate(results_methods$RMSE, list(results_methods$Name), mean)[,2]), linetype="dashed", col="grey");

ggsave(file="~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/plots/jplot_RMSE.pdf",
       jplot_RMSE, height=7, width=7.5);

mean_RMSE_ranking <- aggregate(results_methods$RMSE, list(results_methods$Name), mean);
mean_RMSE_ranking <- mean_RMSE_ranking[order(mean_RMSE_ranking[,2]),];

jplot_MAE <- ggplot(data=results_methods, aes(x=Name, y=MAE, group=1)) + 
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean, geom="crossbar", size=0.5, color="grey") +
  geom_jitter(aes(col=Cell), position = position_jitter(height = 0, width = 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = 90, hjust = 1),
        axis.title=element_text(size=14,face="bold")) +
  scale_colour_brewer(type = "qual", palette = 'Dark2') + 
  xlab('Cell-type deconvolution strategy') + 
  geom_hline(yintercept=min(aggregate(results_methods$MAE, list(results_methods$Name), mean)[,2]), linetype="dashed", col="grey");

ggsave(file="~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/plots/jplot_MAE.pdf",
       jplot_MAE, height=7, width=7.5);

mean_MAE_ranking <- aggregate(results_methods$MAE, list(results_methods$Name), mean);
mean_MAE_ranking <- mean_MAE_ranking[order(mean_MAE_ranking[,2]),];

jplot_R2 <- ggplot(data=results_methods, aes(x=Name, y=R2, group=1)) + 
  stat_summary(fun.y=mean, fun.ymin = mean, fun.ymax = mean, geom="crossbar", size=0.5, color="grey") +
  geom_jitter(aes(col=Cell), position = position_jitter(height = 0, width = 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_text(angle = 90, hjust = 1),
        axis.title=element_text(size=14,face="bold")) +
  scale_colour_brewer(type = "qual", palette = 'Dark2') + 
  xlab('Cell-type deconvolution strategy') + 
  geom_hline(yintercept=max(aggregate(results_methods$R2, list(results_methods$Name), mean)[,2]), linetype="dashed", col="grey");

ggsave(file="~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/plots/jplot_R2.pdf",
       jplot_R2, height=7, width=7.5);

mean_R2_ranking <- aggregate(results_methods$R2, list(results_methods$Name), mean);
mean_R2_ranking <- mean_R2_ranking[order(mean_R2_ranking[,2], decreasing = T),];

methods_to_scatterplot <- as.character(unique(all_predictions$Name));
s_palette <- brewer.pal(6,"Dark2")[order(colnames(real_cc_gs)[-1])];

for(m in methods_to_scatterplot){
  
  sub_method <- all_predictions[all_predictions$Name==m,];
  real_df_ordered <- real_cc_gs[match(sub_method$Sample, real_cc_gs$Sample),];
  i <- 1;
  m_plots <- list();
  
  for(c in colnames(real_df_ordered)[-1]){
    
    data_to_plot <- data.frame(real=as.numeric(unlist(real_df_ordered[c])), predicted=as.numeric(unlist(sub_method[c])));
    m_plots[[i]] <- ggplot(data=data_to_plot, aes(x=real,y=predicted)) + geom_point(col=s_palette[i], size=3) + 
    geom_abline(slope=1, intercept=0, linetype=3) +
    theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlab("Real cellular fractions (%)") + ylab("Predicted cellular fractions (%)") +
    labs(title = paste0("Cell type: ", c), subtitle = paste0('Strategy: ', m));
    i <- i+1;
  }
  
  exp <- do.call("arrangeGrob", c(m_plots, nrow=3, ncol=2)); 
  ggsave(file=paste0('~/Desktop/methylation_clock/polycomb_hypothesis/epigenetic_syndromes/cell_composition/plots/', m,'_scatterplots.pdf'), 
         exp, height=15, width=10);
  
}


# According to RMSE and MAE idol_NFB_houseman is the optimal strategy. 

############################################################
################## End of the script  ######################
############################################################

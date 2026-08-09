# ------------------------------------------------------------------------------
#
#  SES output compilation  
#
# ------------------------------------------------------------------------------

# Author:

# Created: 04/28/2026

# Description :


# House Keeping ----------------------------------------------------------------
rm(list=ls())


# Summary stats and plots  -----------------------------------------------------

# The remainder of code is used to create figures 4 and s6.1-2 and summary 
# stats used in the multidimensional beta change results section


# Summary Stats ----------------------------------------------------------------

# Percent change in observed mean pairwise beta
mpw_pct_list <-lapply(c("mod_B","his_B"), function(i){
  time_name = ifelse(i=="mod_B","Contemporary","Native")
  out_list_t$obs[[i]] %>% 
    summarise(across(everything(),~mean(.x)))
})
do.call(function(x,y) 100*(1-x/y),mpw_pct_list)


# Observed delta LCBD summary
summary(out_list_t$obs$`_d_lcbd`)

# ES delta LCBD summary
summary(out_list_t$empirical_es$`_d_lcbd`)

# Beta barplot (Fig 4a) --------------------------------------------------------

# Extract and summarize pairwise beta diversity
mpw_list <-lapply(c("mod_B","his_B"), function(i){
  time_name = ifelse(i=="mod_B","Contemporary","Native")
  out_list_t$obs[[i]] %>% 
    summarise(across(everything(),~mean(.x))) %>% 
    mutate(time = time_name) %>% 
    tidyr::pivot_longer(-time)
})

# Combine and format for plotting
mpw_df <- bind_rows(mpw_list) %>% 
  mutate(
    dim = sub("_mod","",name),
    dim = sub("_his","",dim),
    facet = substr(name,1,3),
    component = sub(".*B", "", name),
    facet =factor(facet,c("tax","fun","phy")),
    time = factor(time, c("Native", "Contemporary"))
  ) %>% 
  filter(component != "total")

# Plot
pair <- ggplot(mpw_df, aes(x = time, y = value, fill = dim)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~facet)+
  theme_classic(base_size = 25)+
  ylab("")+
  theme(legend.position = "none",
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("#BF4354","#CCA6A7",
                               "#4A7B97","#A6B0BD",
                               "#CEA23F","#DCD0A7"))

# Export
ggsave(
  filename = file.path(fig_dir,"mean_pairwise_beta.png"),
  plot = pair,
  width = 10,
  height = 6,
  dpi = 600
)


# Observed plots (Fig 4c) -----------------------------------------------------
plot_metrics <-  c("_d_lcbd","_d_B")
lapply(plot_metrics,function(i){
  # Subset data
  obs_plot_df <- out_list_t$obs[[i]] %>% 
    tibble::rownames_to_column("COMID") %>% 
    tidyr::pivot_longer(-COMID) %>% 
    mutate(name = sub("_d","",name),
           name = sub("rf","",name))
  
  # Set factor levels
  plot_factor <- c("tax_Btotal","tax_Brepl","tax_Brich",
                   "fun_Btotal","fun_Brepl","fun_Brich",
                   "phy_Btotal","phy_Brepl","phy_Brich")
  obs_plot_df$name <- factor(obs_plot_df$name,plot_factor)
  
  # Plot
  obs_box <- ggplot(
    obs_plot_df, 
    aes(x = name, y = value, fill = name,color = name)
  ) +
    geom_boxplot()+
    stat_summary(
      geom = "crossbar", 
      width = 0.75,         
      fatten = 2,         
      color = "black",        
      fun.data = function(x){
        return(c(y = median(x), ymin = median(x), ymax = median(x))) 
      }
    )+
    geom_hline(yintercept=0, linetype="dashed", color = "red",size = 1)+
    theme_classic(base_size = 25)+
    theme(
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 2.5),
      axis.title.y=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank()) +
    scale_y_continuous(breaks = pretty) +
    scale_fill_manual(values = c("#CEA23F","#CEA23F","#CEA23F",
                                 "#BF4354","#BF4354","#BF4354",
                                 "#4A7B97","#4A7B97","#4A7B97"))+
    scale_color_manual(values = c("#CEA23F","#CEA23F","#CEA23F",
                                  "#BF4354","#BF4354","#BF4354",
                                  "#4A7B97","#4A7B97","#4A7B97"))
  
  # Export
  obs_name <- paste0("obs_",sub("_d_","",i),"_box.png")
  ggsave(
    filename = file.path(fig_dir,obs_name),
    plot = obs_box,
    width = 10,
    height = 6,
    dpi = 600
  )
})


# Effect size plots (FIg 4b,d)  -----------------------------------------------
plot_metrics <-  c("_d_lcbd","_d_B")
lapply(plot_metrics,function(i){
  # Subset and format data
  es_plot_df <- out_list_t$empirical_es[[i]] %>% 
    tibble::rownames_to_column("COMID") %>% 
    tidyr::pivot_longer(-COMID) %>% 
    mutate(name = sub("_d","",name),
           name = sub("rf","",name))
  
  # Set factor levels
  plot_factor <- c("tax_Btotal","tax_Brepl","tax_Brich",
                   "fun_Btotal","fun_Brepl","fun_Brich",
                   "phy_Btotal","phy_Brepl","phy_Brich")
  es_plot_df$name <- factor(es_plot_df$name,plot_factor)
  
  # Remove taxonomic
  es_plot_df <- es_plot_df %>% 
    filter(substr(name,1,1) != "t")
  
  # Plot
  es_box <- ggplot(
    es_plot_df, 
    aes(x = name, y = value, fill = name,color = name)
  )+
    geom_boxplot()+
    stat_summary(
      geom = "crossbar", 
      width = 0.75,         
      fatten = 2,         
      color = "black",        
      fun.data = function(x){
        return(c(y = median(x), ymin = median(x), ymax = median(x)))
      }
    )+
    geom_hline(yintercept=0, linetype="dashed", color = "red",size = 1)+
    theme_classic(base_size = 25)+
    theme(
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 2.5),
      axis.title.y=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank()) +
    scale_y_continuous(breaks = pretty) +
    scale_fill_manual(values = c(
      "#BF4354","#BF4354","#BF4354",
      "#4A7B97","#4A7B97","#4A7B97"))+
    scale_color_manual(values = c(
      "#BF4354","#BF4354","#BF4354",
      "#4A7B97","#4A7B97","#4A7B97"))
  
  #Export
  es_name <- paste0("es_",sub("_d_","",i),"_box.png")
  ggsave(
    filename = file.path(fig_dir,es_name),
    plot = es_box,
    width = 6.6,
    height = 6,
    dpi = 600
  )
})

# Pairwise plots (Fig. s6.1-2)--------------------------------------------------

# Extract delta LCBD data
d_lcbd <-final_list$`_d_lcbd`


# Obs correlation
d_lcbd %>% select(
  tax_Btotal,tax_Brepl,tax_Brich,
  fun_Btotal,fun_Brepl,fun_Brich,
  phy_Btotal,phy_Brepl,phy_Brich) %>% 
  cor()

# ES correlation
d_lcbd %>% select(
  tax_Btotal,tax_Brepl,tax_Brich,
  fun_Btotal_es,fun_Brepl_es,fun_Brich_es,
  phy_Btotal_es,phy_Brepl_es,phy_Brich_es) %>% 
  cor() 

# Pairwise comparison plot data

# obs
plot_dat <- d_lcbd %>%
  select(
    tax_Btotal, tax_Brepl, tax_Brich,
    fun_Btotal, fun_Brepl, fun_Brich,
    phy_Btotal, phy_Brepl, phy_Brich
  )

# es
plot_dat <- d_lcbd %>%
  select(
    tax_Btotal, tax_Brepl, tax_Brich,
    fun_Btotal_es, fun_Brepl_es, fun_Brich_es,
    phy_Btotal_es, phy_Brepl_es, phy_Brich_es
  )

# Open PNG device
png(
  filename = file.path(fig_dir,"dlcbd_pairs_plot_es.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300
)

# Smaller text and points for dense 9x9 matrix
par(
  mar = c(1.2, 1.2, 1.2, 1.2),
  cex = 0.35,         # overall scaling
  pch = 16,
  col = "black"
)

# Pairs plot
pairs(
  plot_dat,
  gap = 0.15,
  cex.labels = 1,
  cex.axis = 1
)

# Close device
dev.off()


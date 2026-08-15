#-------------------------------------------------------------------------------
#
#  Origin-based invadedness data prep
# 
#-------------------------------------------------------------------------------

# Author: 

# Created: 04/28/2026

# Description: Uses raw community data to calculate species richness and 
# community invadedness based on geographic origin of nonnative species. Exports 
# file for use as explanatory variable in RDA.


# House keeping ----------------------------------------------------------------
rm(list = ls())

# Directories
com_dir <- "Diversity Input Data"
pred_dir <- "Analysis_data"

# Load packages
library(dplyr)
library(tidyr)

# Load data
com_df <- readRDS(file.path(com_dir,"raw_community_diversity_input.rds")) %>% 
  ungroup()

# Origin-based invadedness -----------------------------------------------------
inv_df <- com_df %>% 
  
  # Assign scale based native status
  mutate(
    native_status = case_when(
      Native8 == T ~ "native",  # native to that huc 8
      Native8 == F & Native2 == T ~ "prov",  # native to region but not to huc 8
      Native2 ==F & NativeCon ==T ~ "reg",  # native continent but not regions
      NativeCon == F ~ "extr" # not native to continent
    )
  )   %>% 
  
  # Pivot data frame to create columns for the number of each type of species
  group_by(HUC_12,COMID,native_status) %>% 
  summarise(richness = n_distinct(Scientific_Name)) %>%   
  pivot_wider(  
    names_from = native_status, 
    values_from = richness, 
    values_fill = 0
  ) %>% 
  
  # Calculate invadedness as proportion of richness made up by nonnative species
  mutate(
    tot_sp = native + prov + reg + extr,  # calculate proportions of species
    invProv = prov/tot_sp,
    invReg = reg/tot_sp,
    invExtr = extr/tot_sp
  ) %>% 
  ungroup() %>% 
  
  # Designate a column for native taxonomic alpha diversity
  rename(nTaxA = native) %>% 
  ungroup() 

# Export for summary stats script
saveRDS(inv_df,file.path(pred_dir,"origin_invaded.rds"))

# Nonnative species tables -----------------------------------------------------

# tables of species by nonnative origin, and number of sites they occupy as
# that type of nonnative species
nn_table_df <- com_df %>% 
  mutate(
    native_status = case_when(
      Native8 == T ~ "native",  
      Native8 == F & Native2 == T ~ "prov",  
      Native2 ==F & NativeCon ==T ~ "reg",  
      NativeCon == F ~ "extr" 
    )
  ) 

# extra-realm species
sp_summary_exotic <- nn_table_df %>%   
  filter(native_status == "extr") %>% 
  group_by(Scientific_Name) %>% 
  summarise(n_sites = n()) %>% 
  mutate(prop_sites = n_sites/1023)

# Regional
sp_summary_reg <- nn_table_df %>%   
  filter(native_status == "reg") %>% 
  group_by(Scientific_Name) %>% 
  summarise(n_sites = n())%>% 
  mutate(prop_sites = n_sites/1023)

# Provincial
sp_summary_tran <- nn_table_df %>%   
  filter(native_status == "prov") %>% 
  group_by(Scientific_Name) %>% 
  summarise(n_sites = n())%>% 
  mutate(prop_sites = n_sites/1023)




#' Crush PCR. Updated 2022-02-10
#' 
#' Loads an excel file with raw qPCR data and quantifies expression
#' @param file The excel file with the raw data
#' @param groups A data frame containing two columns, "Sample.Name" and "group" for each sample. Can be included as an additional sheet named "groups" in the excel file
#' @param baseline A string to identify the group used for normalization
#' @param housekeeping String containing the label (or part) of the houskeeping gene
#' @param rm.replicates If TRUE removes samples where technical replicates are missing.
#' @param rm.outliers Remove outliers (TRUE or FALSE)
#' @return A Critical_Result object with the plot, summary data and the stats
#' @import ggplot2
#' @import readxl
#' @import stats
#' @export 

crush_PCR <- function(file, 
                    groups = NULL, 
                    baseline = "c", 
                    housekeeping = "GAPDH", 
                    rm.replicates = TRUE, 
                    rm.outliers = FALSE) {
  
  library(readxl)
  library(tidyverse)
  
  dataset <- read_excel(file,sheet=1, range = "A9:J104", 
                        col_names = c("Well", "Sample.Name","Target.Name","Task","Reporter",
                                      "Quencher","RQ","RQmin","RQmax","Ct"))
  
  if(is.null(groups)){
    if("groups" %in% excel_sheets(file)){
      groups<-read_excel(file, sheet="groups")
    }
    else(groups<-data.frame(Sample.Name=unique(dataset$Sample.Name), group="c"))
  }
  groups$group<-as.factor(groups$group)
  
  dataset <- dataset %>% 
    filter(!is.na(dataset$Sample.Name)) %>% 
    select("Sample.Name", "Target.Name","Ct") %>%   
    mutate(Ct = as.numeric(Ct)) %>% 
    group_by(Sample.Name, Target.Name) %>% 
    mutate(mean = mean(Ct), sd = sd(Ct)) %>% 
    ungroup() %>% 
    mutate(dif = abs(mean - Ct))
  
  # Defining gene names
  INT <- unique(dataset$Target.Name[dataset$Target.Name != housekeeping] )
  dataset$Target.Name[dataset$Target.Name == housekeeping] <- "HK"
  dataset$Target.Name[dataset$Target.Name != "HK"] <- "INT"
  
  # Managing samples with lost technical replicates
  if(rm.replicates == TRUE){
    
    faulty <- dataset %>% 
      filter(is.na(Ct)) %>% 
      select(Sample.Name) %>% unique() %>% as.vector()
    
    if(length(faulty$Sample.Name) != 0){
      cat(paste("Samples:", faulty, "\n Have missing technical replicates. \n Removing them from the analysis. \n \n \n"))
    } else{
      cat("All samples have all technical replicates. \n \n \n")
    }
  } else {
    dataset <- dataset %>% filter(!is.na(Ct))
  }
  
  # Managing deviation in technical replicates
  dataset <- rbind(
    dataset %>% 
      filter(sd > 0.5) %>%
      group_by(Sample.Name, Target.Name) %>% 
      filter(dif != max(dif)) %>% 
      mutate(mean = mean(Ct), sd = sd(Ct)) %>% 
      ungroup(),
    dataset %>% 
      filter(sd < 0.5)
  )
  
  # Calculate delta Cts
  dataset <- dataset %>% 
    group_by(Sample.Name, Target.Name) %>% 
    summarise(mean = mean(Ct), sd = sd(Ct)) %>% 
    pivot_wider(
      names_from = Target.Name,
      values_from = c("mean", "sd")) %>% 
    merge(groups, by = "Sample.Name") %>% 
    mutate(dCt = mean_INT - mean_HK,
           ddCt = dCt - mean(dCt[group == baseline]),
           two.ddCt = 2^-ddCt) 
  
  # Removing outliers
  if (rm.outliers == TRUE){
    dataset <- dataset %>% 
      group_by(group) %>% 
      mutate(Q3 = quantile(two.ddCt)[4], Q1 = quantile(two.ddCt)[2], intq = Q3 - Q1) %>% 
      mutate(lim.up = Q3 + intq, lim.dw = Q1 - intq) %>% 
      filter(lim.dw < two.ddCt & two.ddCt < lim.up) %>% 
      select(1:9)
  } else {
    datset <- dataset
  }
  
  # Plotting
  plot <- ggplot(dataset, aes(x = group, y = two.ddCt))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.1)+
    theme(aspect.ratio = 1.62)+
    ylab("normalized expression")+
    xlab("")+
    labs(title = INT)
  
  return(list(dataset, plot))
}


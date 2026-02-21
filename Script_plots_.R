---
  title: "Script plots"
output: html_document

---
  
  ```{r}
library(Seurat)
library(Matrix)
library(future)
library(future.apply)
library(magrittr)
library(SeuratData)
library(patchwork)
#library(batchelor)
library(SingleCellExperiment)
#library(scran)
#library(scMerge)
#library(BiocSingular)
library(ggplot2)
```

```{r}


# save plots
output_path <- "E:/Review/grafici + stat/grafici/"

```

```{r}
#Objects INPUT 
# - obj integrated
# - Annotation file for object4

int_obj_integ_v5 <- readRDS("E:/Review/outputs/outputsint_obj_integ_v5_quarto.rds") 

# annotation for object4
annotazione_obj4 <- read.csv("E:/Review/grafici + stat/GSE161947_iETX_cell_type_annotation (2).csv") 

```
```{r}


# creation of a disctionarry for correct name of the plots
change_names <- c(
  "2022_Jo" = "Amadei et al. 2022 scRNA seq",
  "JH_old_2022" = "Tarazi et al. 2022",
  "JH_new_formatted" = "Yilmaz et al. 2025",
  "Jo_formatted_2021" = "Amadei et al. 2021",
  "obj5_CX_2022" = "Amadei et al. 2022 inDrop seq"
)

```

# ANNOTATION
## Adding annotation for object4

```{r, annotaizone_mancante}
# Adding annotation missing from merged object

# creating the column I need to add the data 
colnames(annotazione_obj4)[1] <- "cells" 
annotazione_obj4$rownames <- paste0( "Jo_formatted_2021_", annotazione_obj4$cells , "-" , annotazione_obj4$type ) # ricreo uguali a quelli dell'oggetto seurat 
rownames(annotazione_obj4) <- annotazione_obj4$rownames


# creating a column to work on
int_obj_integ_v5@meta.data$tmp <- "tmp"

# Adding cell names
int_obj_integ_v5@meta.data$tmp[ grep( "Jo_formatted_2021_", rownames(int_obj_integ_v5@meta.data))] <- sub("[abc]$", "", rownames(int_obj_integ_v5@meta.data)[ grep( "Jo_formatted_2021_", rownames(int_obj_integ_v5@meta.data))] )

# Adding annotation for obj4
int_obj_integ_v5$cell_annotation[ int_obj_integ_v5$cell_annotation == "obj4" ] <- annotazione_obj4[ int_obj_integ_v5@meta.data$tmp[ int_obj_integ_v5$cell_annotation == "obj4"], ]$annot

# Adding "type" column
int_obj_integ_v5$type[ int_obj_integ_v5$group == "Jo_formatted_2021" ] <- annotazione_obj4[ int_obj_integ_v5@meta.data$tmp[ int_obj_integ_v5$group == "Jo_formatted_2021"], ]$type

#check
is.na(int_obj_integ_v5$cell_annotation) %>%  sum() # top non ci sono più buchi vuoti

```

## Uniform macro annotation

```{r, annotazione_3_macro}

# copying cell annotation in macro annotation 
int_obj_integ_v5$macro_annotation <- int_obj_integ_v5$cell_annotation


# manual modifcations of some ambiguous names
int_obj_integ_v5$macro_annotation[ int_obj_integ_v5$macro_annotation %in% c( "Definitive endoderm", "Allantois/ExE-Mesoderm", "ExE-Mesoderm", "ExE VE") ] <- "new_Epiblasto"


# assigning VE tag
int_obj_integ_v5$macro_annotation[ grep( "Visceral|Parietal|endoderm|VE|XEN ETX|Yolk-sac|ExE VE", int_obj_integ_v5@meta.data$cell_annotation) ] <- "new_VE"

# assigning ExE tag
int_obj_integ_v5$macro_annotation[ grep( "Extra-Embryonic ectoderm|Extraembryonic ectoderm|ExE|Chorion|Endocrine|TS compartment|Ectoplacental cone", int_obj_integ_v5@meta.data$cell_annotation) ] <- "new_ExE"

# assigning epiblast tag
int_obj_integ_v5$macro_annotation[! int_obj_integ_v5$macro_annotation %in% c("new_VE",  "new_ExE") ] <- "new_Epiblasto"

# check 
table( int_obj_integ_v5$macro_annotation )

# Ultimate check
int_obj_integ_v5$cell_annotation[ int_obj_integ_v5$macro_annotation == "new_ExE" ] %>% table()
int_obj_integ_v5$cell_annotation[ int_obj_integ_v5$macro_annotation == "new_VE" ] %>% table()
int_obj_integ_v5$cell_annotation[ int_obj_integ_v5$macro_annotation == "new_Epiblasto" ] %>% table()
```


```{r, fig.width=8}

# fig.width=8 fixing plot's width
# fig.high fixing plot's height

# check for the interested annotations

# int_obj_integ_v5$cell_annotation %>% table() # calling all cell's annotations

annotazione_da_plottare <- c( "Yolk-sac/Hindgut", "Amnion", "Yolk-sac/Endothelial", "XEN ETX") # calling only the cell's we're interested in


DimPlot( subset(int_obj_integ_v5, subset = cell_annotation %in% annotazione_da_plottare ),
         reduction = "umap",
         group.by = "cell_annotation"  ) +
  labs(title = "Titolo del plot") # plot title
```

```{r, fig.width=8}

# plot to check the cell's mapping from the datasets

# int_obj_integ_v5$group %>% table() 

dataset_da_plottare <- c( "2022_Jo", "JH_new_formatted", "JH_old_2022","Jo_formatted_2021", "obj5_CX_2022" ) # dataset to plot


DimPlot( subset(int_obj_integ_v5, subset = group %in% dataset_da_plottare ),
         reduction = "umap",
         #cols = c( "#80B1D3","#FAB95B", "#9E3B3B","#B3DE69","#FB8072") , # change colors
         group.by = c( "group", "cell_annotation", "macro_annotation")[1] ) + 
  labs(title = "Titolo del plot")
```


```{r, fig.width=8}


colors <-c( "#8DD3C7","#FFFFB3","#B7A3E3","#FB8072","#80B1D3", "#FDB462","#B3DE69", "#9E3B3B","#D9D9D9","#BC80BD","#CCEBC5","#E9B63B","#576A8F","#FBB4AE","#696FC7", "#75B06F", "#000000" )



DimPlot(int_obj_integ_v5, 
        reduction = "umap",  
        group.by = c("macro_annotation", "seurat_clusters", "macro_annotation")[2], 
        cols = colors ) 

```



## Natural and Synth 

```{r, separazione Nat e synth}

# creating a new metadata
int_obj_integ_v5@meta.data$nat_syn <- "none"

#int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$group == "2022_Jo"] <- int_obj_integ_v5$orig.ident[ int_obj_integ_v5@meta.data$group == "2022_Jo" ]
int_obj_integ_v5@meta.data$nat_syn[  int_obj_integ_v5@meta.data$group == "2022_Jo"] <- int_obj_integ_v5$type[  int_obj_integ_v5@meta.data$group == "2022_Jo"]

int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$group == "JH_old_2022"] <- int_obj_integ_v5$Nat_Syn[ int_obj_integ_v5@meta.data$group == "JH_old_2022" ]
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$group == "JH_new_formatted"] <- int_obj_integ_v5$type[ int_obj_integ_v5@meta.data$group == "JH_new_formatted" ]
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$group == "Jo_formatted_2021"] <- int_obj_integ_v5$type[ int_obj_integ_v5@meta.data$group == "Jo_formatted_2021" ]
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$group == "obj5_CX_2022"] <- int_obj_integ_v5$sample_resource[ int_obj_integ_v5@meta.data$group == "obj5_CX_2022" ]

int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "synthetic" ] <- "stembryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "natural" ] <- "embryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "Synthetic" ] <- "stembryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "Natural" ] <- "embryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "NE" ] <- "embryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "ETiX" ] <- "stembryo"

int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "iETX_4" ] <- "stembryo" 
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "E6_5" ] <- "embryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "E5_5" ] <- "embryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "ETX_4" ] <- "stembryo"

int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "SEM" ] <- "stembryo"
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "TFSEM" ] <- "stembryo"

int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "ExUt" ] <- "embryo" 
int_obj_integ_v5@meta.data$nat_syn[ int_obj_integ_v5@meta.data$nat_syn == "InUt" ] <- "embryo"


int_obj_integ_v5@meta.data$nat_syn %>%  table()
is.na(int_obj_integ_v5@meta.data$nat_syn) %>% sum() # nessun buco 
table(int_obj_integ_v5$group, int_obj_integ_v5$nat_syn)
```

# PLOTS
## Defining colors for datasets and genes



```{r}
# UMAPs espresione dei geni NAT vs Synt per dataset 

genes <- c("Dnmt1", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Uhrf1" )


# per chiarezza visiva associo un colore fisso ad un database (modificabile):
color_dict <- c(
  "2022_Jo" = "#8F0177",
  "JH_old_2022" = "#FDB462",
  "JH_new_formatted" = "#31694E",
  "Jo_formatted_2021" = "#9E3B3B",
  "obj5_CX_2022" = "#4A70A9"
)
```


## Dotplots - gene's expression per dataset 
### v1 - cycle on datasets 

```{r}

save <- TRUE 

dataset_da_plottare <- c( "2022_Jo", "JH_old_2022", "JH_new_formatted", "Jo_formatted_2021", "obj5_CX_2022")  #selecting dataset



# adding metadata
int_obj_integ_v5$macro_annotation_nat_syn <- paste(int_obj_integ_v5$macro_annotation, int_obj_integ_v5$nat_syn, sep = "_")


for ( dataset in dataset_da_plottare){ 
  dotplot <- DotPlot(
    object = subset(int_obj_integ_v5, subset = group == dataset ),
    features = genes,
    cols = c("lightgrey", color_dict[dataset]), 
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    group.by = "macro_annotation_nat_syn",
    #split.by = "nat_syn", 
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius",
    scale.min = NA,
    scale.max = NA
  ) + 
    coord_flip() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )  + 
    labs(title = paste0( "Dotplot dataset: ", change_names[ dataset] ),  
         size = "Percent Expressed"                                      
    ) +                                                             
    guides(                                                          
      size = guide_legend(order = 1)                                 
    )
  
  print(dotplot)
  
  
  
  if(save){  
    
    name_png <- paste0( output_path, "DotPlot_", dataset, ".png" ) 
    
    ggsave(dotplot, filename = name_png, bg = "white",
           width = 7, 
           height = 6)   
    
  }
  
} 

```

### v2 - single dotplot
```{r, fig.width=8}


dataset <- c( "2022_Jo", "JH_old_2022", "JH_new_formatted", "Jo_formatted_2021", "obj5_CX_2022")[3] #choose with number


single_dot_plot <- DotPlot(object = subset(int_obj_integ_v5, subset = group == dataset ),
                           features = genes,
                           cols = c("lightgrey", color_dict[dataset]), 
                           col.min = -2.5,
                           col.max = 2.5,
                           dot.min = 0,
                           dot.scale = 6,
                           group.by = "macro_annotation_nat_syn",
                           #split.by = "nat_syn",
                           cluster.idents = FALSE,
                           scale = TRUE,
                           scale.by = "radius",
                           scale.min = NA,
                           scale.max = NA
) + 
  coord_flip() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )  + 
  labs(title = paste0( "Dotplot dataset: " , dataset))


print(single_dot_plot)



# saving image

if(FALSE){  #  solo se TRUE salva le immagini 
  
  name_png <- paste0( output_path, "DotPlot_", dataset, ".png" ) 
  
  ggsave(single_dot_plot, filename = name_png, bg = "white",
         width = 7, # da modificare
         height = 6)  # da modificare 
  
}

```

### v3 - merged dotpolots

```{r}


dotplot_all_datasets  <- DotPlot(  object = int_obj_integ_v5,
                                   features = genes,
                                   cols = c("lightgrey", "#8C00FF"), 
                                   col.min = -2.5,
                                   col.max = 2.5,
                                   dot.min = 0,
                                   dot.scale = 6,
                                   group.by = "macro_annotation_nat_syn",
                                   #split.by = "nat_syn",
                                   cluster.idents = FALSE,
                                   scale = TRUE,
                                   scale.by = "radius",
                                   scale.min = NA,
                                   scale.max = NA
) + 
  coord_flip() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )  + 
  labs(title = paste0( "Dotplot dataset: " , "tutti assieme"))


print(dotplot_all_datasets)


```

# Timepoints separation
## - adding metadata

```{r}

# adding column for timepoint metadata
dataset_time_point <- c( "2022_Jo", "obj5_CX_2022")  # tengo solo i dataset divisibili per time point 


# crating the metadata t_point 
int_obj_integ_v5@meta.data$t_point <- NULL 
int_obj_integ_v5@meta.data$t_point[ int_obj_integ_v5$group %in% dataset_time_point[1]] <- int_obj_integ_v5@meta.data$timepoint[ int_obj_integ_v5$group %in% dataset_time_point[1]] %>% sub("-$", "", .)
int_obj_integ_v5@meta.data$t_point[ int_obj_integ_v5$group %in% dataset_time_point[2]] <- int_obj_integ_v5@meta.data$stage[ int_obj_integ_v5$group %in% dataset_time_point[2]]


# defining dictionary that matches the comparisons 

tpoint_dict <- list( "2022_Jo" = list( '1' = c("N_E65", "S_ietx5"),
                                       '2' = c("N_E75", "S_ietx6"),
                                       '3' = c("N_E85", "S_ietx8")),
                     
                     "obj5_CX_2022" =  list('1' = c("Day6", "E7.5"),
                                            '2' = c("Day8", "E8", "E8.5", "E8.75"))  
)

```

## - Dotplots generation 

```{r}


save <- TRUE 


for ( dataset in dataset_time_point ){ 
  
  for ( i in 1:length(tpoint_dict[[dataset]]) ){ 
    
    dotplot <- DotPlot(
      object = subset(int_obj_integ_v5, subset = group == dataset & t_point %in% tpoint_dict[[dataset]][[i]]),
      features = genes,
      cols = c("lightgrey", color_dict[dataset]), 
      col.min = -2.5,
      col.max = 2.5,
      dot.min = 0,
      dot.scale = 6,
      group.by = "macro_annotation_nat_syn",
      #split.by = "nat_syn", 
      cluster.idents = FALSE,
      scale = TRUE,
      scale.by = "radius",
      scale.min = NA,
      scale.max = NA
    ) + 
      coord_flip() + 
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )  + 
      labs(title = paste0( "Dotplot dataset: ",  change_names[ dataset]), 
           # subtitle = paste0( paste( tpoint_dict[[dataset]][[i]], collapse = "-" ) ), 
           size = "Percent Expressed"                                      
      ) +                                                            
      guides(                                                         
        size = guide_legend(order = 1)                                
      )
    
    print(dotplot)
    
    
    if(save){  
      
      name_png <- paste0( output_path, "DotPlot_", dataset, "_",
                          paste( tpoint_dict[[dataset]][[i]], collapse = "-" ), ".png" )  
      
      ggsave(dotplot, filename = name_png, bg = "white",
             width = 7, 
             height = 6)  
      
    }
    
  } 
} 

```


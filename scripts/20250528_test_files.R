library(flowCore)
library(flowWorkspace)
library(ncdfFlow)
library(tidyverse)

(gs_files <- list.files("input/",
                        pattern = "fcs",
                        full.names = TRUE))

fs <- read.ncdfFlowSet(files = gs_files)
gs <- GatingSet(fs)

params_ref_markers <- c("LD","CD8","CFSE","GPR56","CD39")

markers <- as_tibble(pData(parameters(gs_cyto_data(gs)[[1]])),
          rownames = "id") %>% 
  mutate(id = paste0(id,"S"),
         type = case_when(
           desc %in% params_ref_markers ~ "reference",
           is.na(desc) ~ "instrument",
           TRUE ~ "functional"
         ),
         desc_orig = desc,
         desc = ifelse(is.na(desc),name,desc),
         isTransformed  = name %in% names(gh_get_transformations(gs[[1]])))

.mats <- gs_pop_get_data(gs,"root")
.mats <- fsApply(.mats,exprs,simplify = FALSE)
.mats <- lapply(.mats,function(mat) {
  as_tibble(mat)%>% 
    select(with(markers,name[type!="instrument"])) %>% 
    rename_with(~ with(markers,desc[type!="instrument"]))
})

.mats <- bind_rows(.mats,
                   .id = "name")


library(ggcyto)
ggcyto(gs,aes(x = CD4,
              y = HLADR),
       subset="root")+
  geom_hex()

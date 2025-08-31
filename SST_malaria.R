# ## 2. SOM_function 
# Applies Self-Organizing Maps to cluster correlation maps into different groups.
# 
# **Input:**
# - correlation_map_total: Flattened array of correlation maps (maps × grid points)
# - num_rows_list, num_columns_list: SOM grid configuration (e.g., 1×3, 2×2, 1×5)
# 
# **Output:**
# - SOM results array (correlation maps × clusters)

library(RNetCDF)
require(kohonen)
library(ncdf4)
library(stats)

# workdir <- "path of the datasets"

workdir <-"/home/mengxinp/projects/def-mengxinp/mengxinp/Malaria/open_source_code/"
num_rows_list <- c(3,2,5,3,7)
num_columns_list <- c(1,2,1,2,1)

repeat_num=5 # you can repeat the algorithm several times to test the robustness of the result

correlation_map_pool_name <- "correlation_map_SSTA_vs_malaria_anomaly_eachepiweek_lead_0to52_week_reorganized_forSOM_Rmorethan02"
ncfname <- paste(correlation_map_pool_name, ".nc", sep = "")
ncin <- open.nc(paste(workdir,"/",ncfname, sep = ""))
print.nc(ncin)

variable="__xarray_dataarray_variable__"
correlation_map_pool_total <- t(var.get.nc(ncin,variable))

cluster_list <- array(0, dim=c(dim(correlation_map_pool_total)[1],length(num_rows_list),repeat_num))

for (seed_index in 1:repeat_num){
  
  pass_rate_list_large <- array(0, dim=c(300,length(num_rows_list)))
  modified_pass_rate_list_large <- array(0, dim=c(300,length(num_rows_list)))
  
  for (index_try in 1:length(num_rows_list)){
    set.seed(seed_index*1000)
    num_rows <- num_rows_list[index_try]
    num_columns <- num_columns_list[index_try]
    num_cluster=num_rows*num_columns
    print(num_cluster)  
    ar.supersom <- supersom(scale(correlation_map_pool_total), somgrid(num_rows, num_columns, "hexagonal") ) 
    summary(ar.supersom)
    clusters1 <- ar.supersom$unit.classif
    cluster_list[,index_try,seed_index]=clusters1
  }
}

save_file_name=paste("SOM_result_cluster_list_",correlation_map_pool_name,".nc", sep = "")
save_file_name_with_path <- paste(workdir,"/",save_file_name, sep = "")
nc_cluster_list <- create.nc(save_file_name_with_path)

dim.def.nc(nc_cluster_list, "epiweek_leadtime", dim(correlation_map_pool_total)[1])
dim.def.nc(nc_cluster_list, "num_cluster", length(num_rows_list))
dim.def.nc(nc_cluster_list, "seed", repeat_num)

##  Create three variables, one as a coordinate variable

var.def.nc(nc_cluster_list, "epiweek_leadtime", "NC_INT", "epiweek_leadtime")
var.def.nc(nc_cluster_list, "num_cluster", "NC_DOUBLE","num_cluster")
var.def.nc(nc_cluster_list, "seed", "NC_DOUBLE","seed")
var.def.nc(nc_cluster_list, "cluster_list","NC_INT",c("epiweek_leadtime","num_cluster","seed"))

##  Define variable values
myepiweek_leadtime        <- c(1:dim(correlation_map_pool_total)[1])
mynum_cluster <- c((num_rows_list[1]*num_columns_list[1]):(num_rows_list[length(num_rows_list)]*num_columns_list[length(num_columns_list)]))
myseed        <- c(1:repeat_num)

# dim(mynum_cluster) <- c(5,2)

##  Put all of the data:  
var.put.nc(nc_cluster_list, "epiweek_leadtime", myepiweek_leadtime)
var.put.nc(nc_cluster_list, "num_cluster", mynum_cluster)
var.put.nc(nc_cluster_list, "seed", myseed)
var.put.nc(nc_cluster_list, "cluster_list", cluster_list)
print(paste(correlation_map_pool_name,"  finished",sep = "  "))










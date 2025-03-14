# SST_malaria_paper
The source code and exampled datasets for the under-reviewed paper entitled "Beyond ENSO: Harnessing Tropical Ocean Variability in Long-lead Vector-borne Disease Prediction" by Pan et al.

We provides two main functions for analyzing relationships between sea surface temperature anomalies (SSTA) and disease incidence, and identify the dynamic SST index:

# 1. Dataset
- SSTA field data: the sea surface temperature anomaly fields in 1st epiweek of 1999 to 17th epiweek of 2023
    - For each grid, the SST anomaly (SSTA) is calculated by subtracting the seasonally varying climatology from the raw SST value to remove the seasonal cycle. 
    - To match the period of the malaria dataset, we reorganized the daily SSTA field into each epidemiological week. 
    - For each grid, the linear trend of the SSTA in each calendar epidemiological week is removed. 
- Disease incidence anomaly time series: malaria anomaly time series over entire loreto region during 1st epiweek of 2000 to 14th epiweek of 2023 (the up-to-date malaria data we have)
    - the malaria anomaly is used, as we focus on the year-to-year variation. It is calculated by subtracting the malaria annual cycle 
    - (i.e., mean malaria incidence in each calendar epidemiological week during 2000-2022) from the total malaria incidence.

Note: To improve the correlation robustness, we concatenate the data for every four epidemiological weeks in the correlation analysis, so you can see the "add" dimension in both SST and Malaria data. It concatenates the data in four sequencial epiweeks. 

# 2. Correlation_map_function (python)
Generates correlation maps between SSTA fields and disease incidence (e.g., malaria) anomaly time series.

**Input:**
- SSTA field data
- Disease incidence anomaly time series

**Output:**
- Correlation map array containing r-values, p-values, slopes, and intercepts

with correlation_map_function, you can generate correlation maps for SSTA and disease incidences in different calendar weeks and with different time lags.
for example, in our study 2756 correlation maps are generated. we group them together and form the "correlation_map_total", as the input of the second function "SOM_function" function.
Multiple correlation maps can be combined into a `correlation_map_total` array for subsequent clustering analysis.
correlation maps with longitude and latitude need to be flattened into one dimension (the grids dimension in SOM_function)

# 2. Self-organzing map function (R)
Applies Self-Organizing Maps to cluster correlation maps into different groups.

**Input:**
- correlation_map_total: Flattened array of correlation maps (maps × grid points)
- num_rows_list, num_columns_list: SOM grid configuration (e.g., 1×3, 2×2, 1×5)

**Output:**
- SOM results array (correlation maps × clusters)

After identifying different clusters of correlatioin maps, the SST monitoring regions are defined as the largest continuous regions with high Malaria-SST correlation. As a result, a dynamic SST index is identified as the remote predictor for malaria. The “dynamic” implies the SST monitoring regions change among different SOM clusters. When applying to other vector-borne diseases, the temporal continuity should also be checked before identifying the SST monitoring regions.

# Model Card for Self-Organizing Map (SOM)

## Model Details
- **Algorithm**: Self-Organizing Map (Kohonen network)
- **Implementation**: `kohonen` R package (v3.0.12)
- **Hyperparameters**:  
  - Grid: 2×2 hexagonal topology
## Intended Use
- **Purpose**: Cluster correlation maps between SSTA and malaria anomalies.  
- **Domain**: Climate-health modeling in tropical regions.  

## Training Data
- **Input**: 2756 correlation maps (flattened to 1D vectors).  
- **Preprocessing**: Min-max normalized to [0, 1].  
- **Source**: Processed SSTA and malaria data (Loreto, Peru).  

## Limitations
- Sensitivity to grid size; tested configurations (1×3, 2×2, 1×5).  
- Assumes stationarity in SSTA-malaria relationships.  

## Ethical Considerations
- **Data Bias**: Malaria incidence data may underrepresent remote areas.  
- **Generalizability**: Validated for Loreto; applicability to other regions untested.  

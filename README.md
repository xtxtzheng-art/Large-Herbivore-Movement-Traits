# Broad-scale patterns of large herbivore species richness are linked to movement-related traits
This repository contains the R scripts and data structures for analyzing the movement traits of large herbivores. The workflow integrates multi-scale environmental data with movement traits to investigate their influence on species richness through a series of statistical steps, including VIF selection, Spatial COM-Poisson modeling, and Independent Explanatory Power (IEP) calculation.



## 📌 Project Overview
The analysis is divided into four main functional steps:
1. 
**VIF Selection**: Multi-stage multicollinearity screening while protecting key ecological variables (NPP, PSN, TSN).
2. 
**Model Fitting**: Fitting Spatial COM-Poisson models using the `spaMM` package to account for spatial autocorrelation.
3. 
**IEP Calculation**: Computing the Independent Explanatory Power by comparing Pseudo- values between full and reduced models.
4. 
**Visualization**: Generating high-quality boxplots to show the distribution of explanatory power across multiple spatial scales.



## 📂 Data Structure
To ensure the script runs correctly, please organize your local directory to match the structure provided in the **[`Processed_Data`](https://github.com/xtxtzheng-art/Large-Herbivore-Movement-Traits.git)** folder of this repository.
### Key Components:
* **Root Level Files**:
* 
`SRMoveTraits_CC.csv`: Contains Species Richness (SR) target variables for different movement traits.
* 
`macro_CC.csv`: Contains geographic coordinates (, ) and fixed macro-environmental variables.
* 
`processed_macro_CC_*.csv`: Macro-scale data specifically processed for various spatial grains (e.g., 10km, 40km, 400km).
* **Trait-Specific Folders**:
Each subfolder (e.g., `dist_homerange`, `maxspeed`) represents a movement trait configuration and must contain:
* 
`movementtraits_id_CC.csv`: The movement trait data for that specific group.
* 
`processed_lm_*_local.csv`: Local-scale environmental predictors.
* 
`processed_lm_*_context*.csv`: Landscape context predictors calculated at matching scales.
* **VIF Input Directory**:
* A separate directory containing the raw CSV files for initial multicollinearity screening (Step 1).

> **Note**: The script uses dynamic path joining. As long as your folder names and file names match the structure in `Processed_Data`, you only need to update the `project_root` variable in the script.



## 🛠 Prerequisites

### R Packages
Ensure you have the following packages installed:
* **Modeling**: `spaMM`, `car`, `MASS`, `spdep`
* **Data Processing**: `dplyr`, `parallel`, `reshape2`
* **Visualization**: `ggplot2`

### Hardware Note
The Spatial COM-Poisson modeling (`corrHLfit`) is computationally intensive and memory-demanding. It is recommended to run the parallel processing on a machine with at least **32GB RAM** and a multi-core CPU (e.g., 18 cores used in the script).



## 🚀 How to Use
1. **Clone the Repository**:
```bash
git clone https://github.com/xtxtzheng-art/Large-Herbivore-Movement-Traits.git

```
2. **Configure Paths**:
Open the R script and modify the `project_root` variable at the top of the file to point to your local data folder.
3. **Run the Analysis**:
Execute the script sequentially. You can adjust the `n_reps` parameter in the `main_loop` to perform a smaller number of iterations for testing purposes.



## 📝 Citation
If you use this code or data in your research, please cite:
* **Study Title**: Broad-scale patterns of large herbivore species richness are linked to movement-related traits
* **Author**: Xueting Zheng/Nanjing University
* **Journal**: Journal of Animal Ecology, 2026

---

## 📧 Contact
For questions regarding the data or code, please open an issue in this repository or contact **xuetingzheng@smail.nju.edu.cn**.


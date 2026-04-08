#  TMRM: A Temporal Multimodal Deep Learning Framework for 24-Hour Activity Cycle Modeling and Major Cardiometabolic Disease Risk Prediction

We developed a unified deep learning framework termed **TMRM (Temporal Multimodal Risk Model)** to model the **24-hour activity cycle (HAC)**—including sedentary behavior (SB), light physical activity (LPA), moderate-to-vigorous physical activity (MVPA), and sleep—and predict the risk of **major cardiometabolic diseases (MCD)**, including obesity, hypertension, dyslipidemia, and type 2 diabetes (T2D).

Unlike conventional approaches that focus primarily on total activity volume, TMRM explicitly captures the **temporal structure and behavioral composition of the full 24-hour cycle**, enabling a more precise characterization of lifestyle-related disease risk.

This framework is embedded within a broader analytical pipeline that integrates **epidemiological modeling, temporal representation learning, causal inference, and behavioral optimization**, providing a complete and modular solution for:

- High-resolution temporal behavior representation across the 24-HAC  
- Multimodal disease risk prediction using deep learning (TMRM) 
- Modeling nonlinear and dose–response relationships  
- Mediation and causal inference analyses  
- Personalized and population-level activity optimization strategies
  
- ---

## Requirements

The main requirements are listed below:

### Python ≥ 3.8
- NumPy  
- Pandas  
- Scikit-learn  
- PyTorch  
- Joblib  
- SciPy  

### R ≥ 4.0
- data.table  
- ggplot2  
- survival  
- splines  
- mediation  
- TwoSampleMR  
- patchwork  
- cowplot  
- scales  

---

## Description of the Source Code

### 1. `Temporal trajectory.R`
Constructs **24-hour behavioral trajectories** across populations.

Includes:
- Time-resolved activity profiling    
- Visualization of temporal patterns  

---

### 2. `TMRM.ipynb`
Implements the **Temporal Multimodal Risk Model (TMRM)**, a deep learning framework combining:

- Dual-branch LSTM (weekday vs weekend)
- Self-attention mechanisms

Used for:
- MCD risk prediction  
- Model evaluation
- Individual-level risk estimation  

---

### 3. `Mediation analysis.R`
Performs **mediation analysis** to quantify indirect effects of behaviors via intermediate biomarkers.

Includes:
- Natural indirect effects  
- Bootstrapped confidence intervals  
- Pathway decomposition  

---

### 4. `24-HAC_cox.ipynb`
Performs large-scale **Cox proportional hazards regression** to evaluate associations between:

- 24-hour activity features (192-dim representation)
- Incident MCD 

Includes:
- Parallel computation  
- Multiple testing correction (FDR)  
- Feature-level hazard ratio estimation  

---

### 5. `Dose-response relationships.R`
Analyzes **nonlinear dose–response relationships** between activity behaviors and disease risk.

Key components:
- Restricted cubic splines (RCS)  
- Threshold and plateau detection  
- Stratified analyses (e.g., sex, age)  

---

### 6. `Mendelian randomization analysis.R`
Conducts **Mendelian Randomization (MR)** to assess causal effects of activity behaviors on disease risk.

Methods:
- IVW (Inverse Variance Weighted)  
- MR-Egger  
- Sensitivity analyses  

---

### 7. `Optimized 24-HAC profile.R`
Generates **population-level and individual-level optimized activity profiles**.

Core functionality:
- Risk-stratified behavioral targets  
- Personalized recommendation algorithm  
- Constraint-based normalization (time-use feasibility)  
- Comparison of observed vs optimized patterns  

---

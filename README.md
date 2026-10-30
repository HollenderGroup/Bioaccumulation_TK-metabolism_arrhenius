# Bioaccumulation_TK-metabolism_arrhenius

This model integrates **time-varying chemical exposure** with **phase I biotransformation** into a **one-compartment toxicokinetic (TK) model**. A validation script for **fluctuating temperatures** (e.g., mesocosm studies) using the **Arrhenius equation** is also provided.

---

## 📄 Associated Publication

This repository accompanies the journal article:

> Raths et al. (2026). *From laboratory to mesocosm: Temperature and species dependent toxicokinetics of the fluorinated herbicide Diflufenican.* DOI: _TBA_

---

## 🔍 Overview

The core of this repository is an **R implementation** of a one-compartment TK model that includes:

- Phase I biotransformation (TK-metabolism)
- Dynamic (time-varying) exposure
- A validation script for temperature correction via the Arrhenius equation

The model was originally developed in MATLAB using the [BYOM platform](https://www.debtox.info/byom.html).

---

## 📁 Repository Contents

| File | Description |
|------|-------------|
| `TK_metabolism_arrhenius_calibration.R` | Model calibration (TK-metabolism with Arrhenius temperature correction) |
| `TK_metabolism_validation_mesocosm.R` | Model validation under fluctuating temperatures (e.g., mesocosm data) |
| `Arrhenius_results.xlsx` | Arrhenius fit results (input for validation, if required) |

---

## 🚀 Getting Started

### 1. Prepare your input data

Make sure you have all required input data:

- Tissue concentrations
- Water concentrations
- Water temperatures

### 2. Calibrate the model

Open `TK_metabolism_arrhenius_calibration.R` and follow the instructions in the script.

> ⚠️ **Important:** Don't forget to update the **molar weights** of your analytes!

### 3. Validate the model

Open `TK_metabolism_validation_mesocosm.R` and follow the instructions in the script.

> ℹ️ An input file containing the Arrhenius fits (`Arrhenius_results.xlsx`) may be needed.

---

## 📚 Citation

If you use this model in your work, please cite:

```
Raths et al. (2026). From laboratory to mesocosm: Temperature and species dependent
toxicokinetics of the fluorinated herbicide Diflufenican. DOI: TBA
```

---

## 🙏 Acknowledgements

Based on the one-compartment TK-metabolism model developed within the [BYOM platform](https://www.debtox.info/byom.html) (MATLAB).

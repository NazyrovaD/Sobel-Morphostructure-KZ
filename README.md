# Sobel-Morphostructure-KZ

**Automated morphostructural analysis using Sobel operator and DEM data for East Kazakhstan**

This repository contains the Google Earth Engine (GEE) scripts and auxiliary materials used for automated detection of morphostructural lineaments based on the Sobel gradient operator.

---

## 🔍 Description
The workflow implements the Sobel operator for regional morphostructural zoning using DEMs (MERIT, SRTM, and ALOS).  
It includes:
- Gradient computation (Sobel, Laplacian, and Canny filters)
- Thresholding of top 15% gradient values
- Orientation and azimuthal analysis
- Enrichment analysis relative to ore deposits
- Visualization and export steps

---

## 📁 Repository Structure
/scripts/ → Google Earth Engine (GEE) JavaScript files
  - `SOBEL.js` — DEM-based Sobel gradient computation and orientation analysis  
  - `ENRICHMENT.js` — Quantitative enrichment and distance analysis of deposits within Sobel zones
/figures/ → Output figures and gradient maps
/data/ → Reference fault and deposit datasets (metadata only)
README.md → This documentation file
LICENSE → MIT License file

---

## 🧭 Study Reference
This repository supports the publication:  
> *Nazyrova D., Kayukov P., Berkinbaev G., Askarov S., Temirbekov N., Temirbekov A., Kasenov S., Temirbekova L.*  
> **Automated Sobel-Based Morphostructural Analysis and Metallogenic Interpretation of the East Kazakhstan Orogenic Province** (Submitted to *MDPI Geosciences*, 2025).

---

## 🧠 Citation
If you use this workflow, please cite the article above or acknowledge this repository.

---

## 📫 Contact
**Corresponding author:**  
Dilyara Nazyrova  
📧 *nazyrovadilyara@gmail.com*  
ECOSERVICE-S LLP, Kazakhstan



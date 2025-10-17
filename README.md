# Sobel-Morphostructure-KZ

**Automated morphostructural analysis using Sobel operator and DEM data for East Kazakhstan**

This repository provides the Google Earth Engine (GEE) scripts and supporting materials used in the morphostructural analysis of East Kazakhstan based on the Sobel gradient operator.

---

## 🧭 Description

The workflow implements the Sobel operator for regional morphostructural zoning using DEMs (MERIT, SRTM, and ALOS).

It includes:
- Gradient computation (Sobel, Laplacian, and Canny filters)
- Thresholding of top 15% gradient values
- Orientation and azimuthal analysis
- Validation against known fault structures
- Enrichment analysis relative to ore deposits
- Visualization and export stages

---

## 🗂️ Repository Structure

├── scripts/
│ ├── SOBEL.js # Main Sobel morphostructural workflow
│ ├── ENRICHMENT.js # Enrichment analysis of deposits vs Sobel zones
│
├── figures/
│ ├── Figure1_AOI_Map.png # Study area (East Kazakhstan)
│ ├── Figure4a_Sobel_Gradient.png # Sobel gradient magnitude
│ ├── Figure4b_Sobel_TopZone.png # Top-15% Sobel gradient mask
│ ├── Figure4c_Sobel_Orientation.png # Orientation map (0–180°)
│ ├── Figure4d_Orientation_Rose.png # Rose diagram of Sobel lineaments
│ ├── Figure5_Enrichment_Curve.png # Enrichment curve (Deposits vs Sobel zones)
│
├── LICENSE
└── README.md


---

## 🖼️ Figures

| File | Description |
|------|--------------|
| `Figure1_AOI_Map.png` | Study area and regional context of East Kazakhstan. |
| `Figure4a_Sobel_Gradient.png` | Sobel gradient magnitude derived from MERIT DEM (90 m). |
| `Figure4b_Sobel_TopZone.png` | Binary mask of the top 15% gradient magnitudes (active morphostructural zones). |
| `Figure4c_Sobel_Orientation.png` | Orientation (0–180°) of Sobel-derived lineaments. |
| `Figure4d_Orientation_Rose.png` | Rose diagram summarizing azimuthal distribution of lineaments. |
| `Figure5_Enrichment_Curve.png` | Enrichment curve showing relationship between deposits and Sobel zones. |

---

## 🧩 How to Cite

If you use this code or materials in your research, please cite:

> “Automated Sobel-Based Morphostructural Analysis and Metallogenic Interpretation Using DEM Data for East Kazakhstan”, *MDPI Geosciences*, 2025 (in review).

---

## ⚙️ License

This repository is released under the **MIT License**.  
You are free to use, modify, and distribute it with proper attribution.

---

## 🌐 Contact

**Corresponding author:** nazyrovadilyara@gmail.com




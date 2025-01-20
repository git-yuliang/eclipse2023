# eclipse2023

These Jupyter Notebook documents and scripts demonstrate the processing and analysis of observational data collected during the total solar eclipse of 2023.

# Input data:

The SDO/AIA data, the $H_{\alpha}$ image observed by the Learmonth Global Oscillation Network Group (GONG), the $TiO\ 7057\ \text{Å}$ image obtained using the Educational Adaptive-optics Solar Telescope (EAST) at the Shanghai Astronomy Museum, China, is availiable at \url{https://github.com/git-yuliang/eclipse2023/tree/main/input/EAST/7057}.
The white-light coronal data from the 2023 total solar eclipse, captured using an iPhone, are available on Figshare via the private link: \url{https://figshare.com/s/91053c9b87dcff703be1} and will be made public upon article acceptance.
Given the substantial volume of white-light coronal data, the file paths were not stored in the default directory (\texttt{./input/*}) but were instead saved in the local environment. In the code, we specified the path as \texttt{/Volumes/WD1T/share/TSE2023/iPhone/}. To rerun the results, please download the dataset and update the path to match your local directory structure.

# Manuscript Title:

# **High-Frequency Quasi-Periodic Oscillations in the Solar Corona Observed with High-Frame-Rate Imaging During the 2023 Total Solar Eclipse**

### **Data Processing and Analysis Notebook**

**Authors**: [Yu Liang]
**Affiliation**: [Shanghai Astronomical Observatory, Chinese Academy of Sciences, Shanghai 200030, China]
**Contact**: [yuliang@shao.ac.cn]
**Last Updated**: [2024-12-27]

---

### **Notebook Overview**

The primary objectives of this Notebook include:

1. Preprocessing of raw white-light corona imaging data (dark field subtraction, flat field correction, and temporal registration).
2. Application of wavelet analysis to detect high-frequency quasi-periodic oscillations (HFQPOs).
3. Visualizations of the solar corona and statistical results for Sun-as-a-star analysis.

---

### **Dependencies**

The analysis relies on the following Python libraries (compatible versions):

- `numpy==1.2.43`
- `pandas==1.5.1`
- `scikit-image==0.20.0`
- `matplotlib==3.7.1`
- `astropy==5.1.1`
- `opencv-python==4.7.0.72`
- `Pillow==9.2.0`
- `scipy==1.13.1`

<!-- A complete list of dependencies can be found in the accompanying `requirements.txt`. -->

---

### **Data Source**

- **Observation Site**: Learmonth Airport, Australia
- **Instrument**: iPhone with 565 nm filter, sampling at 240 fps, frams width * height: 720.0 * 1280.0, image shape:(1280, 720, 3)
- **Date**: April 20, 2023

---

### **Main Structure**

- **Section 1**: Data loading and inspection
- **Section 2**: Preprocessing (dark field, flat field, and registration)
- **Section 3**: Wavelet analysis and signal detection
- **Section 4**: Results visualization and discussion

# **Sun-as-a-star analysis of coronal high-frequency quasi-periodic oscillations observed by an iPhone during the total solar eclipse of 20 April 2023 with 240 fps imaging**

### **Data Processing and Analysis Notebook**

**Authors**: [Yu Liang]
**Affiliation**: [Shanghai Astronomical Observatory, Chinese Academy of Sciences, Shanghai 200030, China]
**Contact**: [yuliang@shao.ac.cn]
**Last Updated**: [2024-12-27]

---

### **Notebook Overview**

The primary objectives of this Notebook include:

1. Preprocessing of raw white-light corona imaging data (dark field subtraction, flat field correction, and temporal registration).
2. Application of wavelet analysis to detect high-frequency quasi-periodic oscillations (HFQPOs).
3. Visualizations of the solar corona and statistical results for Sun-as-a-star analysis.

---

### **Dependencies**

The analysis relies on the following Python libraries (compatible versions):

- `numpy==1.2.43`
- `pandas==1.5.1`
- `scikit-image==0.20.0`
- `matplotlib==3.7.1`
- `astropy==5.1.1`
- `opencv-python==4.7.0.72`
- `Pillow==9.2.0`
- `scipy==1.13.1`

<!-- A complete list of dependencies can be found in the accompanying `requirements.txt`. -->

---

### **Data Source**

- **Observation Site**: Learmonth Airport, Australia
- **Instrument**: iPhone with 565 nm filter, sampling at 240 fps, frams width * height: 720.0 * 1280.0, image shape:(1280, 720, 3)
- **Date**: April 20, 2023

---

### **Main Structure**

- **Section 1**: Data loading and inspection
- **Section 2**: Preprocessing (dark field, flat field, and registration)
- **Section 3**: Wavelet analysis and signal detection
- **Section 4**: Results visualization and discussion

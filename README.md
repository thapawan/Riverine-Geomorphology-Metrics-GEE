## Riverine Geomorphology Metrics using Google Earth Engine
This repository contains the code for a comprehensive analysis of riverine geomorphology. The project uses satellite imagery from the Landsat and Sentinel missions to calculate a range of hydraulic, geometric, and environmental metrics for key river segments over a 25-year period (2000–2024).

### 📜 Project Overview
The core objective of this project is to quantify changes in river channels over time, providing valuable data for hydrological research, environmental monitoring, and water resource management. The analysis is performed every five years for both dry and wet seasons to capture seasonal variability.

The script automates the following workflow:

Image Acquisition: Access and harmonize Landsat and Sentinel data from the years 2000, 2005, 2010, 2015, 2020, and 2024.

Water Classification: Employ multiple water indices, including NDWI, MNDWI, and AWEI, to accurately classify water and non-water pixels.

Centerline Extraction: Use the Medial Axis Transform to derive the river's centerline from the classified water masks.

Metric Calculation: Compute a suite of metrics at both the point and line level, including:

Channel Width (m)

Migration Rate (m)

Curvature

Sinuosity

Stream Power (W/m)

Dam Distance (m)

Runoff (mm)

EVI (Enhanced Vegetation Index)

Data Export: Output the calculated metrics into a structured CSV format for further analysis and visualization.

### 🛠️ Technologies Used
Google Earth Engine (GEE): Cloud-based platform for planetary-scale geospatial analysis.

Python: The primary programming language for scripting the workflow.

geemap: A Python package for interactive mapping with GEE.

scikit-image: Used for the Medial Axis Transform.

Pandas: For data manipulation and metric summarization.

### 📝 How to Use the Code
Prerequisites: Ensure you have access to a Google Earth Engine account and have the necessary Python libraries installed.

Authentication: Authenticate your Google Earth Engine account by running ee.Authenticate() in your environment.

Run the Script: Execute the main Python script (main.py) to start the analysis. The script will automatically connect to GEE, process the data for the specified rivers and time periods, and save the results to your Google Drive.

### 🤝 Contribution
Feel free to open issues or submit pull requests to improve the script, add new features, or expand the analysis to other rivers.

📄 License
This project is licensed under the MIT License.

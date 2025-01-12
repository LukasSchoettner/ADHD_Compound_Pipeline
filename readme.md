# ADHD Compound Pipeline

This project analyzes natural compounds and their relationships to ADHD-related data using various computational tools and workflows. The project integrates a database, R scripts, Python scripts, and Nextflow workflows to ensure reproducibility and efficient processing.

---

## **Project Overview**
This pipeline supports:
1. GEO dataset analysis for ADHD-related studies.
2. Integration with natural compound datasets.
3. Automated workflows for data extraction, transformation, and analysis.
4. A graphical user interface (GUI) built with Streamlit for user-friendly interactions.

---

## **Project Structure**
```plaintext
ADHD_Compound_Pipeline/
├── scripts/                # Python and R scripts for data processing
│   ├── gui_app.py          # Streamlit GUI app
│   ├── etl_script.py       # ETL script for database operations
│   ├── analyze_geo.R       # R script for GEO dataset analysis
│   └── ...
├── logs/                   # Log files for debugging
├── results/                # Output results from the pipeline
├── config.yaml             # Centralized configuration file
├── environment.yml         # Conda environment specification
├── nextflow.config         # Nextflow configuration for workflow management
├── README.md               # Project documentation
└── Dockerfile              # Containerization for reproducibility
```

---

## **Setup Instructions**

### Prerequisites
- Miniconda or Anaconda installed on your system.
- PostgreSQL for database management.

### Step 1: Clone the Repository
```bash
git clone <repository_url>
cd ADHD_Compound_Pipeline
```

### Step 2: Set Up Conda Environment
```bash
conda env create -f environment.yml
conda activate ADHD_Compound_Pipeline
```

### Step 3: Configure the Project
Update the `config.yaml` file with your database connection details and pipeline parameters:
```yaml
# Example
geo_db: "postgresql://user:password@localhost/geo_service_db"
```

### Step 4: Set Up the Database
Ensure your PostgreSQL instance is running and the database schemas are created.

### Step 5: Run the GUI
```bash
streamlit run scripts/gui_app.py
```
Access the GUI in your browser at: `http://localhost:8501`

---

## **Key Files**

### 1. `config.yaml`
Contains configurations for database connections, pipeline parameters, and file paths.

### 2. `etl_script.py`
Extracts data from the database, transforms it, and loads results back into a reporting database.

### 3. `gui_app.py`
A Streamlit-based GUI for running the pipeline interactively.

### 4. `nextflow.config`
Manages the workflows for GEO analysis, molecular docking, and pathway enrichment.

### 5. `environment.yml`
Specifies the Python and R dependencies required for the project.

---

## **Common Commands**

### Update the Conda Environment
```bash
conda env update --file environment.yml --prune
```

### Export the Environment
```bash
conda env export > environment.yml
```

### Run the Pipeline (Nextflow)
```bash
nextflow run main.nf
```

---

## **Contributing**
Feel free to submit issues or pull requests to improve the project.

---

## **License**
This project is licensed under the MIT License. See `LICENSE` for details.

---

### Molecular Docking Automation

**Dependencies:**
- Python 3.8+
- Open Babel (`sudo apt install openbabel`)
- AutoDock Vina
- AutoDockTools
- PyMOL (optional, for visualization)

**Usage:**
1. Upload your ligand and protein files through the GUI.
2. Run the molecular docking pipeline.
3. Check results in the `docking_results/` directory.


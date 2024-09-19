# Functional Profiling of Metagenomic Samples

This project processes metagenomic samples by summarizing the functional categories of proteins based on COG (Clusters of Orthologous Groups) annotations. It generates both non-normalized and normalized data outputs for functional profiling.

The script takes sample `.tsv` files, COG definitions, and COG functional categories files as input and outputs two CSV files:
- `non_normalized.csv`: Contains the raw counts of functional categories.
- `normalized.csv`: Contains the relative abundance (normalized counts) of functional categories.

## Requirements

- **Python 3.6+**
- **Pandas** (for data manipulation)

Install the required Python packages using:

```bash
pip install pandas
```

## Directory Structure

```
.
├── data
│   ├── samples
│   │   ├── sample1_protein-cog.tsv
│   │   ├── sample2_protein-cog.tsv
│   │   └── ...
│   ├── cog-24.def.tab
│   ├── cog-24.fun.edited.tab
├── functional_profiling_script.py
└── README.md
```

### Input Files

- **Sample `.tsv` files**: These files contain the mapping between proteins and their associated COGs. Each file should follow this structure:

```
Protein    COG
S3_00015    COG:COG0190
S3_00016    COG:COG0826
S3_00017    COG:COG1538
...
```

- **COG Definitions File (`cog-24.def.tab`)**: A tab-delimited file containing the mapping of COG IDs to functional categories and other metadata. Expected columns:
  - `COG ID`, `COG Functional category ID`, `COG name`, `Gene`, `Pathway`, `PubMed ID`, `PDB ID`

- **COG Functional Categories File (`cog-24.fun.edited.tab`)**: A tab-delimited file mapping functional category IDs to their groups and descriptions. Expected columns:
  - `COG Functional category ID`, `Functional group`, `RGB color`, `FC description`

### Output Files

- **`non_normalized.csv`**: Contains raw counts of functional categories for each sample.
- **`normalized.csv`**: Contains the relative abundance (normalized counts) for each sample.

## Usage

### Command-line Arguments

The script uses the `argparse` library to handle command-line inputs. Here are the options available:

- `-i` or `--input_directory`: Directory containing sample `.tsv` files (required).
- `-o` or `--output_directory`: Directory to save the output CSV files (required).
- `-d` or `--cog_definitions_file`: Path to the COG definitions file (`cog-24.def.tab`) (required).
- `-f` or `--cog_functional_categories_file`: Path to the COG functional categories file (`cog-24.fun.edited.tab`) (required).

### Example Usage

```bash
python functional_profiling_script.py \
  -i /path/to/samples \
  -o /path/to/output \
  -d /path/to/cog-24.def.tab \
  -f /path/to/cog-24.fun.edited.tab
```

In this example:
- `/path/to/samples` contains `.tsv` files for each metagenomic sample.
- `/path/to/output` is the directory where the results (`non_normalized.csv`, `normalized.csv`) will be saved.
- `/path/to/cog-24.def.tab` is the COG definitions file.
- `/path/to/cog-24.fun.edited.tab` is the COG functional categories file.

### Sample Output

After running the script, two CSV files will be generated in the specified output directory:

#### `non_normalized.csv` (raw counts)

```
Functional group,FC description,sample1,sample2
Cellular processes and signaling,Cell wall/membrane/envelope biogenesis,5,3
Cellular processes and signaling,Defense mechanisms,7,8
Metabolism,Amino acid transport and metabolism,6,2
...
```

#### `normalized.csv` (relative abundance)

```
Functional group,FC description,sample1,sample2
Cellular processes and signaling,Cell wall/membrane/envelope biogenesis,0.176,0.150
Cellular processes and signaling,Defense mechanisms,0.245,0.320
Metabolism,Amino acid transport and metabolism,0.090,0.070
...
```

### Example Files

If you'd like to see an example of the expected input files:

#### `sample1_protein-cog.tsv`

```
S3_00015    COG:COG0190
S3_00016    COG:COG0826
S3_00017    COG:COG1538
S3_00019    COG:COG1136
S3_00020    COG:COG0577
...
```

#### `cog-24.def.tab`

```
COG0001    H    Glutamate-1-semialdehyde aminotransferase    HemL    Heme biosynthesis    NaN    2CFB
COG0002    E    N-acetyl-gamma-glutamylphosphate reductase    ArgC    Arginine biosynthesis    NaN    3DR3
COG0003    P    Anion-transporting ATPase, ArsA/GET3 family    ArsA    NaN    NaN    1F48
...
```

#### `cog-24.fun.edited.tab`

```
COG Functional category ID,Functional group,RGB color,FC description
H,3,DCDCFC,Coenzyme transport and metabolism
J,1,FCCCFC,Translation, ribosomal structure and biogenesis
M,2,DCFCAC,Cell wall/membrane/envelope biogenesis
...
```

### Notes

- Ensure the file paths and filenames are correct when running the script.
- The script will create the output directory if it does not exist.
- Use proper tab-delimited `.tsv` format for all input files.

## Contact

For any questions or issues regarding this script, feel free to contact the project maintainer.

---

Let me know if you'd like any adjustments to this README file!
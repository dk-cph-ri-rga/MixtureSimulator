# MixtureSimulator 1.0

MixtureSimulator is a python script to simulate `.fastq.gz` mixtures by combining single-source `.fastq/fastq.gz` profiles based on user-defined mixture parameters.

## Features
- Supports **single-end** and **paired-end** sequencing data.
- Accepts `.fastq` and `.fastq.gz` files.
- Converts `.fastq` to `.fastq.gz` automatically.
- Logs all parameters, environment, and warnings to a **timestamped log file**.
- Validates read-counts and pairing integrity.
- Random subsampling ensures reproducibility of realistic mixtures.

## Installation

Download the script MixtureSimulator1.0.py and place it in `path/to/folder/src/`.

Also install external dependencies, needed to run MixtureSimulator1.0:
```bash
python        version 3.7 or higher
Biopython     version 1.85
```

## Usage
### To run
```bash
path/to/your/folder/src/MixtureSimulator1.0.py -i <input_dir> -o <output_dir> -b <base_line> -r <repetitions> -n <number_ind> -ra <ratios> [--paired]
```
### Arguments
```bash
-h, --help      Show options
--version       Show version

Required:
 -i, --input_dir        Path to the folder containing single-source fastq/fastq.gz files.
 -o, --output_dir       Path to the output folder for writing mixtures and logs.
 -b, --base_line        Number of reads in simulated mixtures.
 -re, --repetitions     Number of repetitions per simulated mixture.
 -n, --number_ind        Number of individuals per mixture.
 -ra, --ratios           Comma-separated list of mixture ratios.

Optional:
--paired                Enable for paired-end data 
```
### Input files
Supported files are `.fastq` and `.fastq.gz` compressed files. Multiple single source files are necessary to create the mixtures. The script will automatically check if enough `.fastq` files are present in the input directory `-i / --input_dir`, to create the mixtures accordingly to the user-specified parameters `-n / --number_ind`.
When the paired-end mode is enabled `--paired`, the names of the `.fastq` files of the forward and reverse reads from each individual must end with `"_R1"` and `"_R2"`. The script will automatically check if every `"_R1"` file has a `"_R2"` pair.
The script will also ensure the files are not empty and have enough reads to generate a mixture according to the user-defined read counts `-b / --base_line`.

### Output files
The compressed `.fastq.gz` files of the simulated mixtures will be gathered in the output directory `-o / --output_dir`, together with a `.log` file containing all the necessary information about the run. The number of mixtures generated will vary depending on the number of individuals present in the input directory `-i / --input_dir` and the number of contributors to the mixture specified by the user `-n / --number_ind`.
If the paired-end mode is enabled `--paired`, each mixture will have two FASTQ compressed files: `"_R1.fastq.gz"` and `"_R2.fastq.gz"`, for forward and reverse reads.

## Author
Developed by:
```bash
BioInformatic Science Support (BISS)
Copenhagen University
Faculty of Health Sciences
Department of Forensic Medicine
Section of Forensic Genetics
http://retmedicin.ku.dk/english/
```

## Support 
For suggestions and comments, please contact Monica Giuffrida (monicagiuffrida@sund.ku.dk).

## License and citation
Licence can be found at: https://github.com/dk-cph-ri-rga/MixtureSimulator/blob/main/LICENSE.md
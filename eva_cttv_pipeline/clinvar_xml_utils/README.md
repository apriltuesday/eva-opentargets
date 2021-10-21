# ClinVar XML utils

## Overview

The **ClinVar XML utils** package contains utilities and classes to parse the [ClinVar XML](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/) (see [README](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/_README)) and convert the records into internal representation via the following classes:
* **ClinVarDataset**. Iterate through records in ClinVar XML dump and convert them into internal ClinVarRecord representation.
* **ClinVarRecord**. Instances of this class hold data on individual ClinVar records (see [README](https://github.com/EBIvariation/eva-opentargets/tree/master/data-exploration/clinvar-variant-types) for an in-depth explanation of ClinVar's data model). 
* **ClinVarTrait**. Represents a single [ClinVar trait](https://github.com/EBIvariation/eva-opentargets/tree/master/data-exploration/clinvar-variant-types#trait-representation) (usually a ``disease``), with the corresponding database and PubMed cross-references.
* **ClinVarRecordMeasure**. This class represents individual [ClinVar record "measures"](https://github.com/EBIvariation/eva-opentargets/tree/master/data-exploration/clinvar-variant-types#variation-representation). Measures are essentially isolated variants, which can be combined into either ``MeasureSets`` (include one or more Measures) or ``GenotypeSets``.
* **ClinVarRecordMeasureHGVS**. used to determine the type of variant (``trinucleotide_repeat_expansion``, ``short_tandem_repeat_expansion`` or ``None``). 

## Installation

The package can be installed running the following command:
````
pip install clinvar_xml_utils
````

#! TBD: If we add additional (optional) additional packages to the installation, add their installation here. For example: "Developing package X... To install xml_clinvar_utils along the tools you need to develop and run tests, run the following command in your virtualenv: 
    ```bash
    $ pip install "clinvar_xml_utils[dev]"
    ```

#! TBD: Execute tests


## Usage

#! TBD: Include documentation on how to use the utils

## Examples

#! TBD: add examples that can be run within the ``data`` folder to explain the workflow of the script. 






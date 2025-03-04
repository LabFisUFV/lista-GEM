# lista-GEM: Genome-scale metabolic reconstruction of _Lipomyces starkeyi_

![GitHub release (with filter)](https://img.shields.io/github/v/release/LabFisUFV/lista-GEM) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8367982.svg)](https://doi.org/10.5281/zenodo.8367982)


## Description

This repository contains the current genome-scale metabolic reconstruction of _Lipomyces starkeyi_ NRRL Y-11557, named **lista-GEM**. The model distributed on this GitHub repository is continuously updated, with the latest releases available [here](https://github.com/LabFisUFV/lista-GEM/releases). To get access to the model associated to the Almeida et al. (2023) publication, use lista-GEM 1.0.0.

## Citation

* If you use lista-GEM please cite the preprint:
  > Almeida et al. _lista-GEM: the genome-scale metabolic reconstruction of Lipomyces starkeyi._ bioRxiv (2023). [DOI: 10.1101/2023.09.25.559328](https://doi.org/10.1101/2023.09.25.559328)

## Keywords

**Utilisation:** experimental data reconstruction; multi-omics integrative analysis;, _in silico_ strain design; model template   
**Field:** metabolic-network reconstruction   
**Type of model:** reconstruction; curated   
**Model source:** [iYali](http://doi.org/10.1038/npjsba.2016.5) and [rhto-GEM](http://doi.org/10.1002/bit.27162)   
**Omic source:** genomics  
**Taxonomic name:** _Lipomyces starkeyi_   
**Taxonomy ID:** [taxonomy:675824](https://identifiers.org/taxonomy:675824)   
**Genome ID:** [insdc.gca:GCA_018804115.1](https://identifiers.org/insdc.gca:GCA_018804115.1)   
**Metabolic system:** general metabolism   
**Strain:** NRRL Y-11557   
**Condition:** aerobic; defined media   


## Installation

You can obtained the model by any of the following methods:
1. If you have a Git client installed on your computer, you can clone the [`main`](https://github.com/LabFisUFV/lista-GEM) branch of the lista-GEM repository.
2. You can directly download [the latest release](https://github.com/LabFisUFV/lista-GEM/releases) as a ZIP file.
3. If you want to contribute to the development of lista-GEM (see [below](#below)), it is best to [fork](https://github.com/LabFisUFV/lista-GEM/fork) the lista-GEM repository to your own Github account.

## Usage

If you want to use the model for your own model simulations, you can use **any software** that accepts SBML L3V1 FBCv3 formatted model files. This includes any of the following:
* MATLAB-based
  * [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN) version 2.8.3 or later (recommended)  
  * [COBRA Toolbox](https://github.com/opencobra/cobratoolbox)

* Python-based
  * [cobrapy](https://github.com/opencobra/cobrapy)  

***RAVEN Toolbox***
```matlab
model = importModel('lista-GEM.xml')
exportModel(model, 'lista-GEM.xml')
```

***cobrapy***
```python
import cobra
model = cobra.io.read_sbml_model('lista-GEM.xml')
cobra.io.write_sbml_model(model, 'lista-GEM.xml')
```

***COBRA Toolbox*** \*
```matlab
model = readCbModel('lista-GEM.xml')
writeCbModel(model, 'lista-GEM.xml')
```
\* note that some annotation might be lost when exporting the model from COBRA Toolbox.

## Contributing

Contributions are always welcome! Please read the [contributing guideline](.github/CONTRIBUTING.md) to get started.


### Contributors

Code contributors are reported automatically by GitHub under [Contributors](https://github.com/LabFisUFV/lista-GEM/graphs/contributors), while other contributions come in as [Issues](https://github.com/LabFisUFV/lista-GEM/issues).


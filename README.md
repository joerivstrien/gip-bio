# GIP

<!-- ![Tests](https://github.com/joerivstrien/compact-bio/actions/workflows/tests.yml/badge.svg) -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7777356.svg)](https://doi.org/10.5281/zenodo.7777356) -->

### Gaussian Interaction Profiler

## Contents

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Python Package](#python-package)
- [License](#licence)
- [Issues](#issues)
- [Citing CompaCt](#citing-compact)


## About

--needs an about description--

For a more complete description of the software and its applications, please refer to the manuscript (see [here](#citing-compact)).

## Installation

-- needs an installation description -- 

**dependencies**
- [Python >= 3.7](https://www.python.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- NEEDS FULL LIST OF REQUIREMENTS

### pip
installation with pip ensures the python-based dependencies are automatically installed alongside CompaCt. However, MCL (and optionally fastrbo) will still need to be manually installed.

    pip install gip-bio

### Installation from repository

    git clone git@github.com:joerivstrien/compact-bio.git
    cd compact-bio
    pip install .

# Usage

## Input Data

-- needs information on input data -- 
- complexome profile
- annotation file

### complexome profiling data

-- description needs to be fixed/completed -- 
A matrix of expression/abundance data, like a complexome profile. These data should be provided as a tab-separated text file, with the first row containing identifiers. An example of a file containing protein protein expression data is available [here](https://github.com/joerivstrien/gip-bio/blob/main/example_data/CRS17_genenames_renamed.tsv)


### protein annotation file

-- description of annotation file should go here -- 
An example mapping file is available [here](https://github.com/joerivstrien/gip-bio/blob/master/example_data/ENTERFILENAMEHERE)


<!-- ### Command Line Arguments
    #################################
    #### commonly useful options ####
    #################################

    -h, --help          show a complete list of all command line arguments
    -v, --version       show program's version

    ##############################################################
    #### more obscure options, ignoring these is usually fine ####
    ##############################################################

    -m [MIN_SEARCH_WEIGHT], --min-search-weight [MIN_SEARCH_WEIGHT]
                        minimum search weight for rank biased overlap, default=0.99
    --th-criterium {percent,best}
                        threshold criterium for determination of reciprocal top hits.
                        default='percent'
    --th-percent [TH_PERCENT]
                        if th-criterium='percent': top % to consider top hits.
                        default=1.
    --include_within    include within-sample interaction scores in clustering input.
                        default=False
    --wbratio [WBRATIO]   ratio of within/between score averages. default=1
    --skip-filter-clusters
                        use if low-scoring clusters should not be filtered
    --annot-ref-th [ANNOT_REF_TH]
                        min fraction of reference group that must be present in
                        cluster. default=0.5
    --annot-mem-th [ANNOT_MEM_TH]
                        min fraction-present of clust members to be counted for
                        annot-ref-th. default=0.25 -->

## Output

-- needs output annotation -- 

## Python Package

Aside from using the CompaCt command line tool, that performs the complete analysis from expression/abundance datasets to annotated species-specific clusters, CompaCt can also be used as a python package. From python CompaCt can be used more flexibly, for example to rerun only specific steps in the analysis with changed parameters. 

## Licence

    GIP -- Gaussian Interaction Profiler
    Copyright (C) 2023 Radboud University Medical Center

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Issues
If you have questions or encounter any problems or bugs, please report them in the [issue channel](https://github.com/joerivstrien/gip-bio/issues).

## Citing GIP

Publication Pending
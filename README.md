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
- [Citing GIP](#citing-gip)


## About

We introduce the Gaussian Interaction Profiler (GIP), a Gaussian mixture modeling-based clustering workflow for complexome profiling data. GIP assigns proteins to a set number of clusters by modeling the migration profile of each cluster. Using bootstrapping, GIP offers a way to prioritize actual interactors over spuriously comigrating proteins.

For a more complete description of the software and its applications, please refer to the manuscript (see [here](#citing-gip)).

## Installation

GIP is implemented as a flexible python package, which requires installing the package and its dependencies from the python package index (pip).

**dependencies**
- [Python 3](https://www.python.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [scikit-learn](https://scikit-learn.org/)
- [scipy](https://scipy.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)

### pip
installation with pip ensures the dependencies are automatically installed alongside GIP.

    pip install gip-bio

### Installation from repository

    git clone git@github.com:joerivstrien/gip-bio.git
    cd gip-bio
    pip install .

# Usage

## Input Data

GIP takes as input a single complexome profiling dataset, consisting of a series of abundance values that represent the fraction of the migration pattern for each detected protein.

- complexome profile
- protein annotation file

### complexome profiling data

A single complexome profiling dataset, consisting of a matrix of expression/abundance data. These data should be provided as pandas DataFrame, with the index row containing protein identifiers. The GIP package contains a function (process_normalise.parse_profile) to load a complexome profile from a tab-separated text (tsv) file. An example of a file containing a complexome profile is available [here](https://github.com/joerivstrien/gip-bio/blob/main/example_data/test_complexome_profile.tsv)


### protein annotation file

To provide additional protein annotations aside from their identifiers to the output tables containing the GIP analysis results, a table can be provided containing these annotations. This table should be provided as a pandas.DataFrame, where the index contains protein identifiers that match those in the provided complexome profile. An example annotation file, in tab-separated (tsv) format is available [here](https://github.com/joerivstrien/gip-bio/blob/master/example_data/human_annotation.tsv)


## Running a complete GIP analysis

Below is an example script to run a complete GIP analysis from start to finish. Please be aware that, when running GIP from a python script it is important to ensure the code is run from within an "if \__name__ == "\__main__":" block, as in the example below. Running the GIP code outside of a block could result in RuntimeErrors related to multiprocessing, as depending on the multiprocessing start-method each time a new process is started by python it reloads the main module.

The complete results of a GIP run are returned as a dictionary from the main function. Which can then be explored, saved etc. Additionally, when using the "clusttable_fn" and "membertable_fn" parameters, the main result tables are saved to files. These contain detailed information of the clusters and clustered proteins, respectively. 

    from gip.main import main
    import gip.process_normalise as prn
    import pandas as pd

    if __name__ == "__main__":
      # parse complexome profile and protein annotation file
      prof = prn.parse_profile('path/to/profile.tsv')
      annot = pd.read_csv('path/to/annot_fn.tsv',sep='\t',index_col=0)

      # set ratio of clusters relative to number of detected proteins
      clust_ratio = 0.5

      # A standard run, using 4 threads for the bootstrapping
      # the main result tables, containing details regarding the clusters and
      # clustered proteins are saved to tsv files.
      # for this example only 4 bootstrap iterations are run to save time.
      gip_results = main(prof, clust_ratio, annot_df=annot, bs_processes=4,
                         n_bootstraps = 4,
                         clusttable_fn='./clusttable.tsv',
                         membertable_fn='./membertable.tsv')

## Output

The main result of a GIP analysis is a set of clusters. The clusters are annotated with a variety of metrics that facilitate easy interpretation and prioritization of clusters likely corresponding to actual protein complexes. An overview of all resulting clusters with these features is provided in the output as a table ('clusttable'), which can optionally be saved to a file, using the "clusttable_fn" parameter.

Similarly, all clustered proteins are also annotated with a number of metrics reflecting their assigned cluster, abundance and the consistency with which they are part of their cluster. All protein members are provided as a separate table ('membertable), which can also optionally be saved to a file using the "membertable_fn" parameter.

For a complete description of the output from a GIP analysis please refer to the documentation of the main function [here](https://github.com/joerivstrien/gip-bio/blob/main/src/gip/main.py)

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
# TimeWarp
# Language: R
# Input: TXT (parameters)
# Output: CSV (optimal parameters)
# Tested with: PluMA 1.0, R 3.2.5

PluMA plugin that takes a set of time series data samples
and performs time-warping, which aligns each sample to follow
the same "internal clock" based on biological events.

The input for the plugin is a TXT file of keyword-value pairs, with
each keyword and its meaning listed below:

training (training data set)
clinical (testing data set)
prefix (datasets to align)

Datasets are in the following format:

(prefix)_(datasetname)_(time).csv

where (time) is the current timepoint.  These consist of multiple lines of the form:

The aligned datasets are then sent to the output CSV file provided to the plugin.

Note the input TSV and CSV files in the example/ directory are not publically available.
A future goal is to make a synthetic data set.  In the meantime however, one may
use it on their own tab-separated input data.

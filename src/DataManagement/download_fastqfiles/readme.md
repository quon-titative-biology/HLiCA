Files downloaded with lattice requires additional prepping (mostly parsing steps) for the metadata to be fed into the cellranger pipeline.

Dataset from Gruen's lab required some renaming to follow the expected cellranger naming format.

# Downloading GoogleSheetMetadata

https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/edit#gid=0

For the format of the google sheet metdata, it should be

From navigation link
https://docs.google.com/spreadsheets/d/<KEY>/edit#gid=<GID>
to
https://docs.google.com/spreadsheets/d/<KEY>/export?gid=<GID>&format=csv

# Additional reformat/renaming for specific datasets
Check the additional folder in download_fastqfiles/ that corresponds with a dataset. Some datasets required renaming of the fastqfiles to fit the 10x naming scheme.

# Create additional csv files for running cellranger

```
prep_samples.py
```

# learnLyncs
A bit of a guide on using [lyncs_io](https://github.com/Lyncs-API/lyncs.io) to read gaugefields into python

# Conda Notes
Use `environment.yml` for linux. Use `mac_environment.yml` for Mac. Use `win_environment.yml` for Windows.

## Install Environment
```conda env create -f environment.yml```
## Activate/Use
```conda activate lyncs```
## Update (w. new packages)

1. Edit `environment.yml`
2. Deactivate conda environment with `conda deactivate`
3. Update conda environment with `conda env update -f=environment.yml`

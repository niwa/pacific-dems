# pacific-dems
Python code and Jupyter notebooks for creating nation-wide DEMs of Pacific countries from LiDAR DEMs where available and FABDEM elsewhere

## Repository structure
* Python scripts are in the `py_scripts` folder
* Jupyter notebooks for creating country DEMs are in the `notebooks_create_dems` folder
* An unversioned data folder is expected for:
   * input data (LiDAR derived DEMs, and FABDEMs)
   * output data (country boundaries, country DEMs)
   * FABDEM data is expected in a `fabdems_folder`
   * Country LiDAR data is expected in a `lidar` folder nested under the `country_name` folder
   * The country dems are created under the country_name` folder in a folder based on the resolution (e.g. `5m_dems` for a 5m resolution).
   
### Data specifics
This is expected to be included in the unversioned data folder in the base folder of this repository. Please update the paths if you want to use a different location. The data folder should countain a `fabdems` folder with all of the unzipped FABDEM tiles covering the country. These can be downloaded at link to FABDEM v1-2. Any country LiDAR data should be included in a folder with the country name.

The FABDEM and any LiDAR data must be downloaded before hand and placed in the expected location. 

### Generating countrywide DEMs
See the notebooks in the `notebok_create_dems` folder for example jupyter notebooks for creating dems at 5m for the specified countrires. The generated DEMs for the PARTneR countires have been uploaded to [NEXUS}(https://nexus.pacificdata.org/#/) and can be access from there. 

### Contribures
Please contact [Rose Pearson](rose.pearson@niwa.co.nz) for any questions about this repository.

### License
This repository is licensed under an MIT licences agreement. 

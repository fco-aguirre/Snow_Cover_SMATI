## This script is for order de data before the main scripts application:


### To apply this script:

**Run in the terminal:**

    python Files_preparation_v2.py -s <source_file> -d <destination_file> -y <year> -n <check_name>
### 1. -s DIR, source_file DIR
  - this is the source directory where your download MODIS file are located 
  - finish this path with the year to be processed, i.e.: /2015
  - also you can put only '.' if you want to use the example files
### 2. -d DIR, destination DIR, 
  - the path to the place to order your files, must be located inside of Data folder in the Script package: '/Users/xxxx/xxxxx/Snow_Cover_SMATI/Data/*<name_project>*/Refelctance_bands'
  - the *<name_project>* is important to keep in order the data structure of the package.
  - so this sub-folder must be created inside the DATA directory before run the script
  - also you can put only 'Brunswick' if you want to use the example files
### 3. -y YEAR, year, the year that you need to process
### 4. -n check_name, checksums name file that has all the record files downloaded

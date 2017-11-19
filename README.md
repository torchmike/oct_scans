
## Usage

1) Get some example .img data to `/data' folder
2) Run `update_A_scan_ROI.m`
3) Run `process_folder_of_OCT_scans.m`
4) This will process all the .img files found from `/data` folder with the associated custom croppings and generate txt files and pngs
5) run bash data/process.sh which will output a tsv of media, mean and standard deviations for each of the samples 

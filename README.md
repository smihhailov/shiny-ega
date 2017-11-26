# shiny-ega
updated error grid analysis for glucose measurement (ega package) implementation in shiny. It includes Parkes (consensus) for type 1 diabetes and type 2 diabetes and Clarkes error grids along with tools to import, analyze and export the data. The package also enables usage of input data in different units.

There are four options for data entry:
1. Test grid: generates 2500 (50x50) datapoints equally separated to form a grid. It is used for demostration purpuses.
2. Example dataset: the dataset provided in ega package.
3. .csv upload: upload data from the user computer. It has a wide range of data entry options.
4. Manual data entry: simple data entry tool based on Ace Editor.

Data analysis includes criteria set by ISO 15197:2013:
1. Number of datapoints
2. Share of datapoints in relevant zones

Data export can be done through either export of the dataset (with zone labels) or export of the graphical representation of the data.

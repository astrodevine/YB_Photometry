# YB
Yellowball Analysis Code
photometry_draw.py is the code that allows the user to draw a polygon over YBs and store photometry results; see the comments at the beginning of the code for the input files and paths that must be updated

photometry_auto.py does the photometry with only user radius as input

compact.py updates the compactness index for a data set that has completed photometry results stored in a csv file, including vertices for coordinates

ExpertPhotom.py is the up to date manual photometry code

Photometry_Flagging.py recreates images from previously conducted photometry and allows users to apply desired flags

Make_Acknowledgements.py makes and formats an acknowledgements list from crowdsourcing submissions and identifies potentially inappropriate submissions for manual review

make_plots.py reads a photometry csv and creates the F12/F8 histograms and color-color KDE plots with crossmatches. Point-plotted color-color plots can be created by uncommenting the end portion of the code

YB_phot_colemanedits_v2.ipynb is the updated student version of the code

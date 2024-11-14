# gutfit

GUTFIT is a program which scan the parameter space of a GUT model in order to fit fermion masses and mixing parameters. 

The gutfit folder contains, type1and2seewsaw_v4_SUSY.py which contains the info about the model used. It contains also experimentalneutrinomassmatrix.py which contains the data parameters used for the fit. 
This version does the scan linearly in all the free parameters, there are other version I can upload where the scan is more generic with a log scan. 

Enter in the main folder /GUTFIT_24/
To run the complete scan run: 
 python3 examples/multinest_v4.py examples/parameter-cards/param_card3sigma_v4.dat -o $OUTPUT_FOLDER_NAME
 
To compute the measure of a certain point in the parameter space run: 
python3 examples/multinest_v4.py examples/parameter-cards/param_card3sigma_v4.dat $OUTPUT_FOLDER_NAME $NAME_OF_PARAMETER_FILE.txt -o $NAMEPLOT.pdf

The benchmark point with which is possible to validate the scan is BP_141124.txt





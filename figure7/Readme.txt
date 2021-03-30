To generate figure a, use the Jupyter notebook at the "heatmap" folder.

To generate figure b, first, get the folder "greedy_sp_data" from the following link:
"https://github.com/yuhaow/sp-intervention/tree/master/perturb_seq/greedy_sp"
and add it to the folders "ges", "greedy_sp", "pc" respectivley.
Then go to each folder and run the scripts "./autosub.sh", "compute_accuracy.R" and "processing_script.R" in order.
Finally, go back to this folder run the script "plot_script.R".

To generate figure c, first copy "greedy_sp_data" folder into the stabs_gsp folder. 
Then run "./autosub.sh" and "report.R" in order, the DAG would be directly printed on the screen. 

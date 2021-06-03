code documenting the analyses of 'Moment by moment tracking of one-shot learning in a naturalistic story reveals the hippocampo-cortical dynamics of memory-based prediction'

the code documents the analyses in the main part of the manuscript, minor details may vary due to data protection and rewriting/renaming for readability

the order of code in the scripts may not correspond to the order in the manuscript, because the paper was restructured in the publication process


In order of the manuscript:
1) Main_Event_Boundaries_Lag.m produces the figures of agreement and does statistics on the advancement of boundary detection (lag)
    - the data is provided as dataBehavioralSegmentation.mat
2) Main_Event_Boundaries_Consensus.m tests for the increase in agreement on the second run
    - the data is provided as dataBehavioralSegmentation.mat

3) The data for the stem-plot is provided in cloze1_{a/b} and cloze2_{a/b} for the prediction experiment and the replication. 
   Fiure1 e is a simple plot of the difference

4) The script Main_Granger_Causality_Of_Gamma_Envelopes.m produces the main contrast of Cortical Predictive Recall (CPR) 
   the result is stored in a structure named: model_eeg. The crucial field is model_eeg.effect that is plotted in Figure 2
   Note: This script depends on the MVGC toolbox (Barnett&Seth)
5) The script Stats_GrangerCausality.m computes the statistics on this effect and does a number of plots for the SI. for instance all electrodes are 
plotted in ROI colors, etc...

6) The function Main_Server_MI2D_Event_Boundary computes the Mutual Information between every channel and the CPR channels.
This takes a long time even for a single subject and run. The window_width and lag can be restricted further, they are too generous!

7) The scipt Stats_MI2D_Event_Boundary_connectivity computes the statistics on the 2D maps by comparing HC to a visual control ROI. There is also phase shuffled data that nicely supports the same clusters, however, those did not make it into the manuscript.
   Note: All MI scripts depend on Robin Ince's Gaussian-Copula Mutual Information tools

8) The scipt Main_Neural_Timecourse_correlation_With_Boundaries.m computes the mean-projected model differences (Figure 4a) and correlates it with the time-course of agreement (Alternatively, the CPR time-course can be locked to event boundaries (Figure4b)

9) The script Main_Neural_Timecourse_Correlation_With_Prediction.m computes the mean-projected model differences (Figure 4a) and correlates it with the change in cloze probability between a sample that has listened to the story and a sample that hasn't.

10) The script Main_Prediction_Locked_MI_connectivity.m computes mutual information between each channel and the CPR channels at different channel offsets

11) The script Stats_Prediction_Locked_MI.m computes the statistics on Hippocampal MI run2 vs. run1,  HC vs. other, HC 2-1 vs. other 2-1, etc. The script also contains the SI analyses of features (either dropped or on their own) 

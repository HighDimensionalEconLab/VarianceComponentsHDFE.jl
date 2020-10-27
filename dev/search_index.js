var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = VarianceComponentsHDFE","category":"page"},{"location":"#VarianceComponentsHDFE","page":"Home","title":"VarianceComponentsHDFE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package runs Kline, Saggio and Sølvsten (KSS) bias correction of variance components in two-way fixed effects models. The link to the repository is this.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The algorithm's purpose is to find two diagonal matrices that store all P_ii and B_ii, and then proceed to compute the bias-correction of the (biased) plug-in estimators. To do so, you will only be required to run an executable file (app) where the user only needs to input the original data (in .csv format) that contains a first identifier (e.g. worker id), a second identifier (e.g firm id) and outcome (e.g. log wage).","category":"page"},{"location":"","page":"Home","title":"Home","text":"With the use of random projections techniques, it is possible to run leave out estimation of variance components on large datasets. To give an idea, for a dataset containing 5 million person year observations, 1.3 million person effects, 90 thousand firm effects the code takes approximately 20 minutes to compute the relevant leave out matrices on a HPC system with 32 cores assigned. For more information on this please refer to Appendix B in the main paper. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores vector of coefficients, the fixed effects for every observation in the leave-out connected set, as well as the diagonal matrices for P_iis and B_iis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [VarianceComponentsHDFE]","category":"page"},{"location":"#VarianceComponentsHDFE.compute_movers-Tuple{Any,Any}","page":"Home","title":"VarianceComponentsHDFE.compute_movers","text":"compute_movers(first_id,second_id)\n\nReturns a vector that indicates whether the first_id' (e.g. worker) is a mover acrosssecondid' (e.g. firms), as well as a vector with the number of periods that each `firstid' appears.\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.get_leave_one_out_set-NTuple{5,Any}","page":"Home","title":"VarianceComponentsHDFE.get_leave_one_out_set","text":"get_leave_one_out_set(y, first_id, second_id, settings, controls)\n\nReturns a tuple with the observation number of the original dataset that belongs to the Leave-out connected set as described in Kline,Saggio, Solvesten. It also provides the corresponding outcome and identifiers in this connected set. At the current version only `controls=nothing' is supported.\n\n\n\n\n\n","category":"method"},{"location":"#About-the-current-version","page":"Home","title":"About the current version","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The bias-correction currently only runs on a model without controls.\nIf the user wants, they can manually preadjust on her own the outcome. For instance, in an AKM context, the user can run first   y = mboxperson effect + mboxfirm effect + Xb + e and feed into the routine y-Xb as the outcome.","category":"page"}]
}

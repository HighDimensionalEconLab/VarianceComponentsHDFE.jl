var documenterSearchIndex = {"docs":
[{"location":"#VarianceComponentsHDFE","page":"Home","title":"VarianceComponentsHDFE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package estimates a two-way fixed effects model and computes its associated variance components using the methodology developed by Kline, Saggio and Sølvsten (KSS). We provide the usual Julia Package as well as an executable/app that can be run in the terminal. For more details about this please see the corresponding Executable section. The link to the repository can be found  here.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please submit any bug or problem here.","category":"page"},{"location":"#About-the-executable/app","page":"Home","title":"About the executable/app","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores the outcome variable, the first and second set of identifiers, and the estimated fixed effects for every observation that belongs to the leave-out connected set as defined in KSS. The algorithm also saves the statistical leverages, Pii, as well the weighting terms defined as Bii in KSS for a given variance component.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#Windows","page":"Home","title":"Windows","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Open up a powershell terminal. We recommend to run powershell as administrator for installation. To do this, open Windows menu, type \"powershell\". Right-click on the powershell, click \"run as administrator\". \nChange the current directory to where you want to install the executable by typing  in the powershell","category":"page"},{"location":"","page":"Home","title":"Home","text":"cd \"desired_installation_path\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Hint : To copy-paste into the terminal use the standard Ctrl+C and paste into the powershell by using right click.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In the powershell, install the latest version by running:","category":"page"},{"location":"","page":"Home","title":"Home","text":"wget https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5/vchdfe--windows-latest.tar.gz -O vchdfe-windows-latest.tar.gz\n\ntar -xvf vchdfe-windows-latest.tar.gz","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note that to be able to use wget on Windows Internet Explorer should have been launched at least once. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Add the installation directory to PATH. This will allow us to run the program everytime without specifying where the program is installed. To do so copy and paste the following line: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"setx PATH \"$env:path;$pwd\\vchdfe\\bin\" -m","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note: This change will be permanent only if you ran powershell as administrator. Otherwise, everytime you need to run the program you need to specify the installation folder : we would have to type  \"installation_path\"\\\\vchdfe\\\\bin\\\\vchdfe instead of vchdfe everytime we want to run the program. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(OPTIONAL) You can test the program using the sample test file provided with the executable:","category":"page"},{"location":"","page":"Home","title":"Home","text":"vchdfe vchdfe\\bin\\test.csv","category":"page"},{"location":"#MacOS","page":"Home","title":"MacOS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Open Terminal: Press COMMAND + SPACE to open spotlight search, and type terminal and hit RETURN.\nTo download the compressed file in some desired installation directory, run the following code:","category":"page"},{"location":"","page":"Home","title":"Home","text":"cd desired_installation_path\nwget -qO- https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.1/vchdfe-v0.1.5.1-macos-latest.tar.gz | tar -xzv","category":"page"},{"location":"","page":"Home","title":"Home","text":"(RECOMMENDED): To add the bin folder to the PATH, you have to modify the .bash_profile file in your home directory. Add the following line to the bottom of .bash_profile and save it. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"export PATH=\"~/vchdfe/bin:$PATH\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"To source the changes in .bash_profile file, type:","category":"page"},{"location":"","page":"Home","title":"Home","text":"source ~/.bash_profile","category":"page"},{"location":"#Executable-guide","page":"Home","title":"Executable guide","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Before this, make sure you've completed the installation steps shown above. To use the executable you only need to open the terminal. Windows users can press the keys Windows + R, then type \"powershell\" and press Enter. MacOS users can press COMMAND + SPACE to open spotlight search, and type terminal and hit RETURN. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The basic syntax of this command is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"vchdfe path_to_data [--option option_value]","category":"page"},{"location":"","page":"Home","title":"Home","text":"where [] denote optional arguments. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note: You can use this simple syntax if you added the binary path to your system path. See the Note in Step 4 of installation steps for more details. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can see a list of examples provided in a section below or type in the powershell","category":"page"},{"location":"","page":"Home","title":"Home","text":"vchdfe --help","category":"page"},{"location":"","page":"Home","title":"Home","text":"to see a complete list of available arguments:","category":"page"},{"location":"","page":"Home","title":"Home","text":"positional arguments:\n  path                  path to CSV file containing data\n\noptional arguments:\n  --first_id FIRST_ID   column index in CSV file for the first ID\n                        The user should specify the most \"granular\"\n                        id here. For instance, in the AKM context,\n                        this id would correspond to worker ids.\n                        (type: Int64, default: 1)\n  --second_id SECOND_ID\n                        column index in CSV file for the second ID\n                        Use the less granular type (e.g. Firm in the \n                        AKM example)\n                        (type: Int64, default: 2)\n  --observation_id OBSERVATION_ID\n                        column index in CSV file for observation (e.g.\n                        Wage). (type: Int64, default: 4)\n  --no_first_id_effects No computing and showing of first_id effects\n  --no_cov_effects      No computing and showing of covariance effects\n  --algorithm ALGORITHM\n                        type of algorithm: Exact or JLA. It defaults\n                        to be Exact if the number of observations is\n                        less than 5000, and JLA otherwise. (default:\n                        \"Default\")\n  --simulations SIMULATIONS\n                        number of simulations in the JLA algorithm. It\n                        defaults to 100 * log(#total fixed effect)\n                        (type: Int64, default: 0)\n  --header              CSV file contains header\n  --first_id_display FIRST_ID_DISPLAY\n                        The display text associated with first_id\n                        (e.g. Person). (default: \"Person\")\n  --second_id_display SECOND_ID_DISPLAY\n                        The display text associated with second_id\n                        (e.g. Firm) (default: \"Firm\")\n  --observation_id_display OBSERVATION_ID_DISPLAY\n                        The display text associated with observable_id\n                        (e.g. Wage) (default: \"Wage\")\n  --detailed_output_path DETAILED_OUTPUT_PATH\n                        path to the CSV for the detailed output for\n                        each observable (default:\n                        \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\variance_components.csv\")\n  --results_path RESULTS_PATH\n                        path to the results of the output (default:\n                        \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\results.txt\")\n  --write_detailed_CSV  write the detailed output to a CSV\n  --write_results       write the results to a file\n  --print_level PRINT_LEVEL\n                        Level of verbosity of output. (type: Int64,\n                        default: 1)\n  -h, --help            show this help message and exit","category":"page"},{"location":"","page":"Home","title":"Home","text":"The executable can provide two different type of outputs. The first (that we refer to as Results) is a log file, in txt format, that writes a summary of the main results : bias-corrected components and some statistics from the leave-out sample. The second type of output (that we refer to as Detailed Output) is a DataFrame, stored in csv format, that saves important objects computed throughout the algorithm such a subset of the original data corresponding to the leave-out sample, the corresponding Pii and Bii, and the vector of coefficients used to compute the plug-in variance components. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The program provides two outputs: a file contains the main results and a detailed output file. ","category":"page"},{"location":"#Results-output","page":"Home","title":"Results output","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"By using the argument --write_results followed by --results_path MYFILE, you can specify that you want to write the results in MYFILE. ","category":"page"},{"location":"#Detailed-output","page":"Home","title":"Detailed output","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can save all the output details in a CSV file. You need to use --write_detailed_CSV followed by --detailed_output_path MYFILE, where MYFILE is the path to the detailed output file. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The detailed output file includes:                                                                                                                                                                                                                              ","category":"page"},{"location":"","page":"Home","title":"Home","text":"observation : observation identifier in the original dataset that belong to the Leave-out connected set. For instance, if the number 3 appears in this vector, it means that the third observation in the original dataset belongs to the leave-out connected set. \nfirst_id: the first identifier corresponding to each observation in obs (e.g. worked ids in the leave-out connected set).\nsecond_id: the second identifier corresponding to each observation in obs (e.g. firm ids in the leave-out connected set).\nbeta: the vector of coefficients used to compute the plug-in variance components.\nD_alpha: the fixed effect for the first identifier corresponding to each observation. \nF_psi: the fixed effect for the second identifier corresponding to each observation. \nPii: statistical leverage corresponding to each observation in obs.\nBii_first: The Bii for the variance of first_id effects corresponding to each observation in obs.\nBii_second: The Bii for the variance of second_id effects corresponding to each observation in obs.\nBii_cov: The Bii for the co-variance of first_id and second_id effects corresponding to each observation in obs.","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Suppose we have a dataset my_data.csv that is stored in \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\\". The structure of the data is such that the fifth column corresponds to the worker identifier and the third column corresponds to the firm identifier. Moreover, we want to store a summary of the results in this location \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\summary.txt\", and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\output.csv\". ","category":"page"},{"location":"","page":"Home","title":"Home","text":"To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the powershell ","category":"page"},{"location":"","page":"Home","title":"Home","text":"vchdfe \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\my_data.csv\" --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\output.csv\" --results_path \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\summary.txt\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"To run the same thing while specifying 1000 simulations for the JLA algorithm that estimates (Pii,Bii) described in the computational appendix of KSS, we type in the powershell ","category":"page"},{"location":"","page":"Home","title":"Home","text":"vchdfe \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\my_data.csv\" --first_id 5 --second_id 3 --write_results  --write_detailed_CSV --detailed_output_path \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\output.csv\" --results_path \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\summary.txt\" --algorithm JLA --simulations 1000","category":"page"},{"location":"","page":"Home","title":"Home","text":"To only obtain the bias-correction for the variance of firm effects, we type in the powershell ","category":"page"},{"location":"","page":"Home","title":"Home","text":"vchdfe \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\my_data.csv\" --first_id 5 --second_id 3 --no_first_effects --no_cov_effects --write_results  --write_detailed_CSV --detailed_output_path \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\output.csv\" --results_path \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\summary.txt\" --algorithm JLA --simulations 1000","category":"page"},{"location":"#Functions-in-this-package","page":"Home","title":"Functions in this package","text":"","category":"section"},{"location":"#Main-Functions","page":"Home","title":"Main Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"    leave_out_estimation(y,first_id,second_id,controls,settings)\n    get_leave_one_out_set(y, first_id, second_id, settings, controls)    ","category":"page"},{"location":"#VarianceComponentsHDFE.leave_out_estimation-NTuple{5,Any}","page":"Home","title":"VarianceComponentsHDFE.leave_out_estimation","text":"leave_out_estimation(y, first_id, second_id, controls, settings)\n\n\nReturns the bias-corrected components, the vector of coefficients, the corresponding fixed effects for every observation, and the diagonal matrices containing the Pii and Biis. \n\nArguments\n\ny: outcome vector\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on VCHDFESettings\ncontrols: at this version only controls=nothing is supported.\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.get_leave_one_out_set-NTuple{5,Any}","page":"Home","title":"VarianceComponentsHDFE.get_leave_one_out_set","text":"get_leave_one_out_set(y, first_id, second_id, settings, controls)\n\n\nReturns a tuple with the observation number of the original dataset that belongs to the Leave-out connected set as described in Kline,Saggio, Solvesten. It also provides the corresponding outcome and identifiers in this connected set. \n\nArguments\n\ny: outcome vector\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on VCHDFESettings\ncontrols: at this version only controls=nothing is supported.\n\n\n\n\n\n","category":"method"},{"location":"#Auxiliary-Functions","page":"Home","title":"Auxiliary Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"    find_connected_set(y, first_idvar, second_idvar, settings)\n    prunning_connected_set(yvec, first_idvar, second_idvar, obs_id, settings)\n    drop_single_obs(yvec, first_idvar, second_idvar,obs_id)\n    compute_movers(first_id,second_id)\n    eff_res(::ExactAlgorithm, X,first_id,second_id,match_id, K, settings)\n    eff_res(lev::JLAAlgorithm, X,first_id,second_id,match_id, K, settings)\n    compute_matchid(second_id,first_id)   ","category":"page"},{"location":"#VarianceComponentsHDFE.find_connected_set-NTuple{4,Any}","page":"Home","title":"VarianceComponentsHDFE.find_connected_set","text":"find_connected_set(y, first_idvar, second_idvar, settings)\n\n\nReturns a tuple of observation belonging to the largest connected set with the corresponding identifiers and outcomes. This requires to have the data sorted by first identifier, and time period (e.g. we sort by worked id and year). This is also the set where we can run AKM models with the original data.\n\nArguments\n\ny: outcome (e.g. log wage)\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.prunning_connected_set-NTuple{5,Any}","page":"Home","title":"VarianceComponentsHDFE.prunning_connected_set","text":"prunning_connected_set(yvec, first_idvar, second_idvar, obs_id, settings)\n\n\nThis function prunes the dataset from articulation points. If the first identifier is worker id it means that it prunes workers that would disconnect the graph if they were dropped.\n\nArguments\n\nyvec: outcome (e.g. log wage)\nfirst_idvar: first identifier (e.g. worker id)\nsecond_idvar: second identifier (e.g. firm id)\nobs_id: observation identifier.\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.drop_single_obs-NTuple{4,Any}","page":"Home","title":"VarianceComponentsHDFE.drop_single_obs","text":"drop_single_obs(yvec, first_idvar, second_idvar, obs_id)\n\n\nThis function drops observations that correspond with first identifiers with a single observation. For example, if first identifier is worker id, it will drop observations for workers that only appear once in the data.\n\nArguments\n\nyvec: outcome (e.g. log wage)\nfirst_idvar: first identifier (e.g. worker id)\nsecond_idvar: second identifier (e.g. firm id)\nobs_id: observation identifier.\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.compute_movers-Tuple{Any,Any}","page":"Home","title":"VarianceComponentsHDFE.compute_movers","text":"compute_movers(first_id, second_id)\n\n\nReturns a vector that indicates whether the first_id (e.g. worker) is a mover across second_id (e.g. firms), as well as a vector with the number of periods that each first_id appears.\n\nArguments\n\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.eff_res-Tuple{ExactAlgorithm,Any,Any,Any,Any,Any,Any}","page":"Home","title":"VarianceComponentsHDFE.eff_res","text":"eff_res(_, X, first_id, second_id, match_id, K, settings)\n\n\nThis function computes the diagonal matrices containing Pii and Bii under the Exact Algorithm. See appendix in KSS for more information.\n\nArguments\n\nX: the design matrix in the linear model.\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nmatch_id: match identifier between first and second identifier.For example, match can be the identifier of every worker-firm combination.\nK: number of covariates in addition to the fixed effects. Currently only 0 is supported.\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.eff_res-Tuple{JLAAlgorithm,Any,Any,Any,Any,Any,Any}","page":"Home","title":"VarianceComponentsHDFE.eff_res","text":"eff_res(lev, X, first_id, second_id, match_id, K, settings)\n\n\nThis function computes the diagonal matrices containing Pii and Bii under Johnson-Linderstrauss Algorithm. See appendix in KSS for more information.\n\nArguments\n\nX: the design matrix in the linear model.\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nmatch_id: match identifier between first and second identifier.For example, match can be the identifier of every worker-firm combination.\nK: number of covariates in addition to the fixed effects. Currently only 0 is supported.\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"#VarianceComponentsHDFE.compute_matchid-Tuple{Any,Any}","page":"Home","title":"VarianceComponentsHDFE.compute_matchid","text":"compute_matchid(second_id, first_id)\n\n\nComputes a match identifier for every combination of first and second identifier. For example, this can be the match identifier of worker-firm combinations.\n\nArguments\n\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\n\n\n\n\n\n","category":"method"},{"location":"#Datatypes-in-this-package","page":"Home","title":"Datatypes in this package","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"    ExactAlgorithm\n    JLAAlgorithm\n    VCHDFESettings\n    ","category":"page"},{"location":"#VarianceComponentsHDFE.ExactAlgorithm","page":"Home","title":"VarianceComponentsHDFE.ExactAlgorithm","text":"struct ExactAlgorithm <: AbstractLeverageAlgorithm\n\nData type to pass to VCHDFESettings type, to indicate Exact algorithm\n\n\n\n\n\n","category":"type"},{"location":"#VarianceComponentsHDFE.JLAAlgorithm","page":"Home","title":"VarianceComponentsHDFE.JLAAlgorithm","text":"struct JLAAlgorithm <: AbstractLeverageAlgorithm\n\nData type to pass to VCHDFESettings type, to indicate JLA algorithm\n\nFields\n\nnum_simulations: number of simulations in estimation. If num_simulations = 0, defaults to 100 * log(#total fixed effect)\"\n\n\n\n\n\n","category":"type"},{"location":"#VarianceComponentsHDFE.VCHDFESettings","page":"Home","title":"VarianceComponentsHDFE.VCHDFESettings","text":"struct VCHDFESettings{LeverageAlgorithm}\n\nThe VCHDFESettings type is to pass information to methods regarding which algorithm to use. \n\nFields\n\ncg_maxiter: maximum number of iterations (default = 300)\nleverage_algorithm: which type of algorithm to use (default = ExactAlgorithm())\nfirst_id_effects: includes first id effects. At this version it is required to include the firstideffects. (default = true)\ncov_effects: includes covariance of first-second id effects. At this version it is required to include the cov_effects. (default = true)\nprint_level: prints the state of the program in std output. If print_level = 0, the app prints nothing in the std output. (default = 1)\nfirst_id_display_small: name of the first id in lower cases (default = person)\nfirst_id_display: name of the first id (default = Person)\nsecond_id_display_small: name of the second id in lower cases (default = firm)\nsecond_id_display: name of the second id (default = Firm)\nobservation_id_display_small: name of the observation id in lower cases (default = wage)\nobservation_id_display: name of the observation id (default = Wage)\n\n\n\n\n\n","category":"type"},{"location":"#About-the-current-version","page":"Home","title":"About the current version","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The bias-correction currently only runs on a two-way fixed effects model without controls.\nThe user can manually preadjust the outcome. For instance, in an AKM context, the user can run first   y = pe + fe + Xb + e , where pe are person effects and fe are firm effects, and feed into the routine y-Xb as the outcome.","category":"page"}]
}

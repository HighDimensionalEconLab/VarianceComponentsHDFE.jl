var documenterSearchIndex = {"docs":
[{"location":"Executable/#Contents","page":"Executable","title":"Contents","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Pages = [\"Executable.md\"]\nDepth = 3","category":"page"},{"location":"Executable/#About-the-executable/app","page":"Executable","title":"About the executable/app","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores the outcome variable, the first and second set of identifiers, and the estimated fixed effects for every observation that belongs to the leave-out connected set as defined in KSS. The algorithm also saves the statistical leverages, Pii, as well the weighting terms defined as Bii in KSS for a given variance component. More details of this are provided below at the Results output and Detailed csv output sections. ","category":"page"},{"location":"Executable/#Installation","page":"Executable","title":"Installation","text":"","category":"section"},{"location":"Executable/#Windows","page":"Executable","title":"Windows","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Download our latest version of the package from the following link. Move this file to desired_installation_path.\nOpen up a powershell terminal. We recommend to run powershell as administrator for installation. To do this, open Windows menu, type \"powershell\". Right-click on the powershell, click \"run as administrator\". \nChange the current directory to to the desired installation path by typing  in the powershell\ncd \"desired_installation_path\"\nHint : To copy-paste into the terminal use the standard Ctrl+C and paste into the powershell by using right click.\nIn the powershell, install the latest version by running:\ntar -xvf vchdfe--windows-latest.tar.gz\n(RECOMMENDED) A folder named vchdfe/bin now exists inside your installation path. Add the path of that folder to PATH. This will allow us to run the program everytime without specifying where the program is installed. \n(OPTIONAL) You can test the program using the sample test file provided with the executable. If you ran the previous step you may run:\nvchdfe vchdfe\\bin\\test.csv\nOtherwise, you will have to run\n.\\vchdfe\\bin\\vchdfe .\\vchdfe\\bin\\test.csv\nor you could test with some other file you store at the bin folder.\nTo set the number of threads used for parallel computing in the code, you need to use the set command before running vchdfe command. This might be very important if you intend to run bias-correction for very large datasets. Before running the program you may set the number of threads to, say 4, by entering the following line in the Powershell:\n$env:JULIA_NUM_THREADS=4\nTypically, you want to set that number to the number of cores in your computer. You can now proceed to close the terminal.","category":"page"},{"location":"Executable/#PATH-in-Windows","page":"Executable","title":"PATH in Windows","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The PATH is an important concept when working on the command line. It's a list of directories that tell your operating system where to look for programs, so that you can just write program instead of some_folders\\program. But different operating systems have different ways to add a new directory to it. ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"You can take a look at your current PATH by running the following line in the Powershell:","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"$env:path","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"We will describe two ways to add our program to PATH.","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The first way consists on running  this line from the installation folder\nsetx PATH \"$env:path;$pwd\\vchdfe\\bin\" -m\nNote: This change will be permanent only if you ran powershell as administrator. You may need to restart the powershell for this action to take place. You can confirm the change by running $env:path again, and checking whether the last part of this string contains vchdfe\\bin.\nIf the previous step didn't work for you, you can still modify the PATH variable manually. To do so you can press (Windows Key + R) to open the Run Dialog, then type\nrundll32.exe sysdm.cpl,EditEnvironmentVariables\nand press the Enter Key. Find the PATH variable, select it, and click Edit. Then you can click New to add a new program path to the PATH list. You will write the path desired_installation_path\\vchdfe\\bin, where desired_installation_path is the folder where you installed the program. ","category":"page"},{"location":"Executable/#MacOS","page":"Executable","title":"MacOS","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Download our latest version of the package from the following link. Move this file to the desired installation path.\nYou can unpack the .tar.gz file automatically by double-clicking the icon. \nA folder named vchdfe/bin now exists inside your installation path. You can double click the icon named vchdfe inside this folder. If the system warns you that it's not a recognized app you must click OK, to give permission.\n(OPTIONAL) You can test the program using the sample test file provided with the executable. You may run the following line:\nvchdfe vchdfe/bin/test.csv\nor you can test with some other file you store at the bin folder. \nTo set the number of threads used for parallel computing in the code, you need to use the set command before running vchdfe command. This might be very important if you intend to run bias-correction for very large datasets. Before running the program you may set the number of threads to, say 4, by entering the following line in the Powershell:\nset JULIA_NUM_THREADS=4\nTypically, you want to set that number to the number of cores in your computer. You can now proceed to close the terminal.","category":"page"},{"location":"Executable/#Syntax-of-the-program","page":"Executable","title":"Syntax of the program","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Before this, make sure you've completed the installation steps shown above. To use the executable you only need to open the terminal. Windows users can press the keys Windows + R, then type \"powershell\" and press Enter. MacOS users can click to open the app manually. You can also see some examples in the Typical Workflow sections below.","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The basic syntax of this command is ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"vchdfe path_to_data [--option option_value]","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"where [] denote optional arguments. ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"You can see a list of examples provided in a section below or type in the powershell","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"vchdfe --help","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"to see a complete list of available arguments:","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"positional arguments:\n  path                  path to CSV file containing data\n\noptional arguments:\n  --first_id FIRST_ID   column index in CSV file for the first ID\n                        The user should specify the most \"granular\"\n                        id here. For instance, in the AKM context,\n                        this id would correspond to worker ids.\n                        (type: Int64, default: 1)\n  --second_id SECOND_ID\n                        column index in CSV file for the second ID\n                        Use the less granular type (e.g. Firm in the \n                        AKM example)\n                        (type: Int64, default: 2)\n  --outcome_id OUTCOME_ID\n                        column index in CSV file for outcome (e.g.\n                        Wage). (type: Int64, default: 4)\n  --no_first_id_effects No computing and showing of first_id effects\n  --no_cov_effects      No computing and showing of covariance effects\n  --algorithm ALGORITHM\n                        type of algorithm: Exact or JLA. It defaults\n                        to be Exact if the number of observations is\n                        less than 5000, and JLA otherwise. (default:\n                        \"Default\")\n  --covariates NAME1 NAME2 ... \n                        Column names of the covariates that will be \n                        partialled out from the outcome. If dataset\n                        contains no headers use Column1, Column2, etc.\n                        to refer to such columns. (default: [])\n  --do_lincom           Perform linear combination inference.\n  --lincom_covariates NAME1 NAME2 ...\n                        Column names of the covariates used in the linear\n                        combination inference step. If dataset\n                        contains no headers use Column1, Column2, etc.\n                        to refer to such columns. (default: [])\n  --simulations SIMULATIONS\n                        number of simulations in the JLA algorithm. If\n                        no number or zero is provided, it will default\n                        to 200 simulations within the algorithm.\n                        (type: Int64, default: 0)\n  --header              CSV file contains header\n  --first_id_display FIRST_ID_DISPLAY\n                        The display text associated with first_id\n                        (e.g. Person). (default: \"Person\")\n  --second_id_display SECOND_ID_DISPLAY\n                        The display text associated with second_id\n                        (e.g. Firm) (default: \"Firm\")\n  --outcome_id_display OUTCOME_ID_DISPLAY\n                        The display text associated with outcome_id\n                        (e.g. Wage) (default: \"Wage\")\n  --detailed_csv_path DETAILED_CSV_PATH\n                        path to the CSV for the detailed output for\n                        each observable (default:\n                        \"current_directory/variance_components.csv\")\n  --results_path RESULTS_PATH\n                        path to the results of the output (default:\n                        \"current_directory/results.txt\")\n  --write_detailed_csv  write the detailed output to a CSV\n  --write_results       write the results to a file\n  --print_level PRINT_LEVEL\n                        Level of verbosity of output. (type: Int64,\n                        default: 1)\n  -h, --help            show this help message and exit","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The executable can provide two different type of outputs. The first (that we refer to as Results) is a log file, in txt format, that writes a summary of the main results : bias-corrected components and some statistics from the leave-out sample. The second type of output (that we refer to as Detailed Output) is a DataFrame, stored in csv format, that saves important objects computed throughout the algorithm such a subset of the original data corresponding to the leave-out sample, the corresponding Pii and Bii, and the fixed effects that can be used to compute the plug-in variance components. ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The program provides two outputs: a file contains the main results and a detailed output file. ","category":"page"},{"location":"Executable/#Results-output","page":"Executable","title":"Results output","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"By using the argument --write_results followed by --results_path MYFILE, you can specify that you want to write the table with the main results in MYFILE. This will store some summary statistics of the leave out connected sample, the bias-corrected variance components. If do_lincom is activated, it will also store the results of the inference part, where it regresses second effects (e.g. firm effects) against some covariates.","category":"page"},{"location":"Executable/#Detailed-csv-output","page":"Executable","title":"Detailed csv output","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"You can save all the output details in a CSV file. You need to use --write_detailed_csv followed by --detailed_csv_path MYFILE, where MYFILE is the path to the detailed output file. ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The detailed output file includes:                                                                                                                                                                                                                              ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"observation : observation identifier in the original dataset that belong to the Leave-out connected set. For instance, if the number 3 appears in this vector, it means that the third observation in the original dataset belongs to the leave-out connected set. \nfirst_id_old: the first identifier of the original dataset corresponding to each observation in obs (e.g. worked ids in the leave-out connected set).\nsecond_id_old: the second identifier of the original dataset corresponding to each observation in observation (e.g. firm ids in the leave-out connected set).\nfirst_id: the first identifier computed in the routine corresponding to each observation in obs (e.g. worked ids in the leave-out connected set). The maximum of this vector corresponds to the number of first effects (e.g. worker fixed effects).\nsecond_id: the second identifier computed in the routine corresponding to each observation in observation (e.g. firm ids in the leave-out connected set). The maximum of this vector corresponds to the number of second effects (e.g. firm fixed effects).\nDalpha: the fixed effect for the first identifier corresponding to each observation. \nFpsi: the fixed effect for the second identifier corresponding to each observation. \nPii: statistical leverage corresponding to each observation in observation. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.\nBii_first: The Bii for the variance of first_id effects corresponding to each observation in observation. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.\nBii_second: The Bii for the variance of second_id effects corresponding to each observation in observation. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.\nBii_cov: The Bii for the co-variance of first_id and second_id effects corresponding to each observation in observation. If the leave out strategy was done at the match level this number is repeated for every observation that belongs to the same match.","category":"page"},{"location":"Executable/#Typical-Executable-Workflow-(Windows)","page":"Executable","title":"Typical Executable Workflow (Windows)","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"You begin by opening the Powershell (as administrator), and typing ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"cd \"path_to_dataset\"","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"At this point we recommend to exploit parallelization of the underlying code by changing an environment variable to be equal to the number of cores in your computer. Let num_cores be the number of cores. You can run the following line","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"$env:JULIA_NUM_THREADS=num_cores","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Then, if you completed all instalation steps you may run","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"vchdfe my_data.csv --OPTIONS","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"where OPTIONS depend on the structure of your data. ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"If you didn't complete step 5 of the installation, you will have to run instead ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"installation_folder\\vchdfe\\bin\\vchdfe my_data.csv --OPTIONS","category":"page"},{"location":"Executable/#Detailed-examples","page":"Executable","title":"Detailed examples","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Suppose we have a dataset my_data.csv that is stored in \"project_path\". The structure of the data is such that the fifth column corresponds to the worker identifier, third column corresponds to the firm identifier, and the fourth column corresponds to the outcome. Moreover, we want to store a summary of the results in this location \"project_path/summary.txt\", and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  \"project_path/output.csv\". ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt\nTo run the same thing while but we want to partial out the covariates in columns 7 and 8, we type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --covariates Column7 Column8--write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt \nTo only obtain the bias-correction and regress the second effects (e.g. firm effects) onto the observables in columns 7 and 8 instead, we type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --do_lincom --lincom_covariates Column7 Column8 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt  ","category":"page"},{"location":"Executable/#Typical-Executable-Workflow-(MacOS)","page":"Executable","title":"Typical Executable Workflow (MacOS)","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"You begin by clicking the app that can be found inside the vchdfe/bin folder. You may need to click OK to give permission to use the app.  ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"At this point we recommend to exploit parallelization of the underlying code by changing an environment variable to be equal to the number of cores in your computer. Let num_cores be the number of cores. You can run the following line","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"set JULIA_NUM_THREADS=num_cores","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Then, if you completed all instalation steps you may run","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"cd path_to_data\nvchdfe my_data.csv --OPTIONS","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"where OPTIONS depend on the structure of your data. ","category":"page"},{"location":"Executable/#Detailed-examples-2","page":"Executable","title":"Detailed examples","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Suppose we have a dataset my_data.csv that is stored in \"project_path\". The structure of the data is such that the fifth column corresponds to the worker identifier, the third column corresponds to the firm identifier and the fourth column corresponds to the outcome. Moreover, we want to store a summary of the results in this location \"project_path/summary.txt\", and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  \"project_path/output.csv\".  You begin by opening the app as shown in the previous section.","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt\nTo run the same thing while but we want to partial out the covariates in columns 7 and 8, we type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --covariates Column7 Column8--write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt \nTo only obtain the bias-correction and regress the second effects (e.g. firm effects) onto the observables in columns 7 and 8 instead, we type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --do_lincom --lincom_covariates Column7 Column8 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt  ","category":"page"},{"location":"Executable/#About-the-current-version","page":"Executable","title":"About the current version","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The bias-correction currently only runs on a two-way fixed effects model without controls.\nThe user can manually preadjust the outcome. For instance, in an AKM context, the user can run first   y = pe + fe + Xb + e , where pe are person effects and fe are firm effects, and feed into the routine y-Xb as the outcome.","category":"page"},{"location":"#VarianceComponentsHDFE","page":"Home","title":"VarianceComponentsHDFE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package estimates a two-way fixed effects model and computes its associated variance components using the methodology developed by Kline, Saggio and Sølvsten (KSS). We provide the usual Julia Package as well as an executable/app that can be run in the terminal. For more details about this please see the corresponding Executable section. The link to the repository can be found  here.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please report any bug or problem here.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"Executable.md\", \"Package.md\"]\nDepth = 3","category":"page"},{"location":"Package/#Contents","page":"Package","title":"Contents","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"Pages = [\"Package.md\"]\nDepth = 3","category":"page"},{"location":"Package/#Setting-up-Development-Enviroment","page":"Package","title":"Setting up Development Enviroment","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"Install Julia and GitHub Desktop - not strictly required but never hurts to have it!\nInstall vscode and follow basic instructions in https://github.com/ubcecon/tutorials/blob/master/vscode.md\nIn particular, https://github.com/ubcecon/tutorials/blob/master/vscode.md#julia, making sure to do the code formatter step.\nand the git settings in https://github.com/ubcecon/tutorials/blob/master/vscode.md#general-packages-and-setup\nClone the repo by either:\nClicking on the Code then Open in GitHub Desktop.\nAlternatively, you can go ] dev https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl in a Julia REPL and it will clone it to the .julia/dev folder.\nIf you wanted to, you could then drag that folder back into github desktop.\nOpen it in vscode by right-clicking on the folder it installed to, and then opening a vscode project.\nOpen the Julia repl in vscode  (Ctrl-Shift-P and then go Julia REPL or something to find it.\ntype ] instantiate to install all of the packages.  Get coffee.\nIn the REPL run ] test and it should do the full unit test.","category":"page"},{"location":"Package/#Functions-in-this-package","page":"Package","title":"Functions in this package","text":"","category":"section"},{"location":"Package/#Main-Function","page":"Package","title":"Main Function","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"    leave_out_KSS(y,first_id,second_id;controls, do_lincom , Z_lincom , lincom_labels , settings)       ","category":"page"},{"location":"Package/#VarianceComponentsHDFE.leave_out_KSS-Tuple{Any,Any,Any}","page":"Package","title":"VarianceComponentsHDFE.leave_out_KSS","text":"leave_out_KSS(y, first_id, second_id; controls, do_lincom, Z_lincom, lincom_labels, settings)\n\n\nReturns a tuple with the observation number of the original dataset that belongs to the Leave-out connected set as described in Kline,Saggio, Solvesten. It also provides the corresponding outcome and identifiers in this connected set. \n\nArguments\n\ny: outcome vector\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\ncontrols: covariates that will be partialled out from outcome before it performs KSS.\ndo_lincom: boolean that indicates whether it runs inference. \nZ_lincom: matrix of covariates to be used in lincom inference.\nlincom_labels: vector of labels of the columns of Z_lincom.\nsettings: settings based on VCHDFESettings\ncontrols: at this version only controls=nothing is supported.\n\n\n\n\n\n","category":"method"},{"location":"Package/#Auxiliary-Functions","page":"Package","title":"Auxiliary Functions","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"    find_connected_set(y, first_idvar, second_idvar, settings)\n    get_leave_one_out_set(yvec, first_idvar, second_idvar, obs_id, settings)\n    leave_out_estimation(y,first_id,second_id,controls,settings)\n    compute_movers(first_id,second_id)\n    compute_matchid(second_id,first_id)  \n    lincom_KSS(y,X, Z, Transform, sigma_i; lincom_labels ) ","category":"page"},{"location":"Package/#VarianceComponentsHDFE.find_connected_set-NTuple{4,Any}","page":"Package","title":"VarianceComponentsHDFE.find_connected_set","text":"find_connected_set(y, first_idvar, second_idvar, settings)\n\n\nReturns a tuple of observation belonging to the largest connected set with the corresponding identifiers and outcomes. This requires to have the data sorted by first identifier, and time period (e.g. we sort by worked id and year). This is also the set where we can run AKM models with the original data.\n\nArguments\n\ny: outcome (e.g. log wage)\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.get_leave_one_out_set-NTuple{5,Any}","page":"Package","title":"VarianceComponentsHDFE.get_leave_one_out_set","text":"get_leave_one_out_set(y, first_id, second_id, settings, controls)\n\n\nReturns a tuple with the observation number of the original dataset that belongs to the Leave-out connected set as described in Kline,Saggio, Solvesten. It also provides the corresponding outcome and identifiers in this connected set. \n\nArguments\n\ny: outcome vector\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on VCHDFESettings\ncontrols: at this version only controls=nothing is supported.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.leave_out_estimation-NTuple{5,Any}","page":"Package","title":"VarianceComponentsHDFE.leave_out_estimation","text":"leave_out_estimation(y, first_id, second_id, controls, settings)\n\n\nReturns the bias-corrected components, the vector of coefficients, the corresponding fixed effects for every observation, and the diagonal matrices containing the Pii and Biis. \n\nArguments\n\ny: outcome vector\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on VCHDFESettings\ncontrols: matrix of control variables. At this version it doesn't work properly for very large datasets.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.compute_movers-Tuple{Any,Any}","page":"Package","title":"VarianceComponentsHDFE.compute_movers","text":"compute_movers(first_id, second_id)\n\n\nReturns a vector that indicates whether the first_id (e.g. worker) is a mover across second_id (e.g. firms), as well as a vector with the number of periods that each first_id appears.\n\nArguments\n\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.compute_matchid-Tuple{Any,Any}","page":"Package","title":"VarianceComponentsHDFE.compute_matchid","text":"compute_matchid(second_id, first_id)\n\n\nComputes a match identifier for every combination of first and second identifier. For example, this can be the match identifier of worker-firm combinations.\n\nArguments\n\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.lincom_KSS-NTuple{5,Any}","page":"Package","title":"VarianceComponentsHDFE.lincom_KSS","text":"lincom_KSS(y, X, Z, Transform, sigma_i; lincom_labels)\n\n\nThis function regresses fixed effects based onto some observables. See appendix in KSS for more information.\n\nArguments\n\ny: outcome variable.\nX: the design matrix in the linear model.\nZ: matrix of observables to use in regression.\nTransform: matrix to compute fixed effects (e.g. Transform = [0 F] recovers second fixed effects).\nsigma_i: estimate of the unbiased variance of observation i.\nlincom_labels: labels of the columns of Z.\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"Package/#Datatypes-in-this-package","page":"Package","title":"Datatypes in this package","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"    ExactAlgorithm\n    JLAAlgorithm\n    VCHDFESettings\n    ","category":"page"},{"location":"Package/#VarianceComponentsHDFE.ExactAlgorithm","page":"Package","title":"VarianceComponentsHDFE.ExactAlgorithm","text":"struct ExactAlgorithm <: AbstractLeverageAlgorithm\n\nData type to pass to VCHDFESettings type, to indicate Exact algorithm\n\n\n\n\n\n","category":"type"},{"location":"Package/#VarianceComponentsHDFE.JLAAlgorithm","page":"Package","title":"VarianceComponentsHDFE.JLAAlgorithm","text":"struct JLAAlgorithm <: AbstractLeverageAlgorithm\n\nData type to pass to VCHDFESettings type, to indicate JLA algorithm\n\nFields\n\nnum_simulations: number of simulations in estimation. If num_simulations = 0, defaults to 100 * log(#total fixed effect)\"\n\n\n\n\n\n","category":"type"},{"location":"Package/#VarianceComponentsHDFE.VCHDFESettings","page":"Package","title":"VarianceComponentsHDFE.VCHDFESettings","text":"struct VCHDFESettings{LeverageAlgorithm}\n\nThe VCHDFESettings type is to pass information to methods regarding which algorithm to use. \n\nFields\n\ncg_maxiter: maximum number of iterations (default = 300)\nleave_out_level: leave-out level (default = match)\nleverage_algorithm: which type of algorithm to use (default = JLAAlgorithm())\nfirst_id_effects: includes first id effects. At this version it is required to include the firstideffects. (default = true)\ncov_effects: includes covariance of first-second id effects. At this version it is required to include the cov_effects. (default = true)\nprint_level: prints the state of the program in std output. If print_level = 0, the app prints nothing in the std output. (default = 1)\nfirst_id_display_small: name of the first id in lower cases (default = person)\nfirst_id_display: name of the first id (default = Person)\nsecond_id_display_small: name of the second id in lower cases (default = firm)\nsecond_id_display: name of the second id (default = Firm)\noutcome_id_display_small: name of the observation id in lower cases (default = wage)\noutcome_id_display: name of the observation id (default = Wage)\n\n\n\n\n\n","category":"type"},{"location":"Package/#Typical-Julia-Workflow","page":"Package","title":"Typical Julia Workflow","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"#Load the required packages\nusing VarianceComponentsHDFE, DataFrames, CSV\n\n#Load dataset\ndata = DataFrame!(CSV.File(\"test.csv\"; header=false))\n\n#Extract vectors of outcome, workerid, firmid\nid = data[:,1]\nfirmid = data[:,2]\ny = data[:,3]\n\n#You can define the settings using our structures\nJL = JLAAlgorithm(num.simulations = 300)\nmysettings = VCHDFESettings(leverage_algorithm = JL, first_id_effects=true, cov_effects=true)\n\n#Run KSS with no controls \nθ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid)\n\n#Create some controls and run the routine where we partial out them\ncontrols = indexin(year,unique(sort(year)))\ncontrols = sparse(collect(1:size(y,1)), controls, 1, size(y,1), maximum(controls))\ncontrols = controls[:,1:end-1]\n\nθ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid; controls)\n\n#Perform Lincom Inference using a Region Dummy\ndata = DataFrame!(CSV.File(\"lincom.csv\"; header=false))\nid = data[:,1]\nfirmid = data[:,2]\ny = data[:,5]\nregion = data[:,4] \nregion[findall(region.==-1)].=0\n\nθ_first, θ_second, θCOV = leave_out_KSS(y,id,firmid; do_lincom = true , Z_lincom = region, lincom_labels = [\"Region Dummmy\"] )\n\n","category":"page"}]
}

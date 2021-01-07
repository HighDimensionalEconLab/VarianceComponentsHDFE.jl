var documenterSearchIndex = {"docs":
[{"location":"Executable/#Contents","page":"Executable","title":"Contents","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Pages = [\"Executable.md\"]\nDepth = 3","category":"page"},{"location":"Executable/#About-the-executable/app","page":"Executable","title":"About the executable/app","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The algorithm prints the plug-in and the bias-corrected variance components estimators for the first identifier effects (e.g. variance of worker effects), the second identifier effects (e.g. variance of firm effects), and the covariance of both these effects (e.g. covariance of worker-firm effects). The user may choose to compute only a subset of these three components. Additionally, the executable will create a CSV file that stores the outcome variable, the first and second set of identifiers, and the estimated fixed effects for every observation that belongs to the leave-out connected set as defined in KSS. The algorithm also saves the statistical leverages, Pii, as well the weighting terms defined as Bii in KSS for a given variance component.","category":"page"},{"location":"Executable/#Installation","page":"Executable","title":"Installation","text":"","category":"section"},{"location":"Executable/#Windows","page":"Executable","title":"Windows","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Download our latest version of the package from the following link. Move this file to desired_installation_path.\nOpen up a powershell terminal. We recommend to run powershell as administrator for installation. To do this, open Windows menu, type \"powershell\". Right-click on the powershell, click \"run as administrator\". \nChange the current directory to where you want to the desired installation path by typing  in the powershell\ncd \"desired_installation_path\"\nHint : To copy-paste into the terminal use the standard Ctrl+C and paste into the powershell by using right click.\nIn the powershell, install the latest version by running:\ntar -xvf vchdfe--windows-latest.tar.gz\n(RECOMMENDED) Add the installation directory to PATH. This will allow us to run the program everytime without specifying where the program is installed. ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"(OPTIONAL) You can test the program using the sample test file provided with the executable. If you ran the previous step you may run:\nvchdfe vchdfe\\bin\\test.csv\nOtherwise, you will have to run\n.\\vchdfe .\\vchdfe\\bin\\test.csv","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"To set the number of threads used for parallel computing in the code, you need to use the set command before running vchdfe command. This might be very important if you intend to run bias-correction for very large datasets. Before running the program you may set the number of threads to, say 4, by entering the following line in the Powershell:\n$env:JULIA_NUM_THREADS=4\nTypically, you want to set that number to the number of cores in your computer. You can now proceed to close the terminal.","category":"page"},{"location":"Executable/#PATH-in-Windows","page":"Executable","title":"PATH in Windows","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The PATH is an important concept when working on the command line. It's a list of directories that tell your operating system where to look for programs, so that you can just write program instead of some_folders\\program. But different operating systems have different ways to add a new directory to it. ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"You can take a look at your current PATH by running the following line in the Powershell:","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"   $env:path\n   ```\n\nWe will describe two ways to add our program to `PATH`.\n\n1. The first way consists on running  this line from the installation folder\n   \n\n   ```\n   setx PATH \"$env:path;$pwd\\vchdfe\\bin\" -m\n   ```\n\n   Note: This change will be permanent only if you ran powershell as administrator. You may need to restart the powershell for this action to take place. You can confirm the change by running `$env:path` again, and checking whether the last part of this string contains `vchdfe\\bin`.\n\n2. If the previous step didn't work for you, you can still modify the `PATH` variable manually. To do so you can press (`Windows Key` + `R`) to open the Run Dialog, then type\n\n   ```\n   rundll32.exe sysdm.cpl,EditEnvironmentVariables\n   ``` \n   and press the `Enter Key`. Find the `PATH` variable, select it, and click `Edit`. Then you can click `New` to add a new program path to the `PATH` list. You will write the path `desired_installation_path\\vchdfe\\bin`, where `desired_installation_path` is the folder where you installed the program. \n\n## MacOS\n\n1. Download our latest version of the package from the following [link](https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl/releases/download/v0.1.5.9/vchdfe-v0.1.5.9-macos-latest.tar.gz). Move this file to the desired installation path.\n2. Open Terminal: Press COMMAND + SPACE to open spotlight search, and type terminal and hit RETURN.\n\n3. You may unpack the .tar.gz file automatically when you double-click the icon. Otherwise, you may run the following code:\n\n    ```\n    cd desired_installation_path\n    ```\n    \n    ```\n    gunzip -c vchdfe-v0.1.5.3-macos-latest.tar.gz | tar xopft -\n    ```\n\n    And then close the terminal.\n\n4.  (RECOMMENDED) Add the installation directory to [PATH](#path-in-macos). This will allow us to run the program everytime without specifying where the program is installed. \n\n\n5. (OPTIONAL) You can test the program using the sample test file provided with the executable. If you ran the previous step you may run:\n\n   ```\n   vchdfe vchdfe/bin/test.csv\n   ```\n\n   Otherwise, you will have to run\n\n   ```\n   ./vchdfe ./vchdfe/bin/test.csv\n   ```   \n6. To set the number of threads used for parallel computing in the code, you need to use the set command before running `vchdfe` command. This might be very important if you intend to run bias-correction for very large datasets. Before running the program you may set the number of threads to, say 4, by entering the following line in the Powershell:\n\n    ```\n    set JULIA_NUM_THREADS=4\n    ```\n    \n   Typically, you want to set that number to the number of cores in your computer. You can now proceed to close the terminal.\n\n\n## PATH in MacOS\n\nThe PATH is an important concept when working on the command line. It's a list of directories that tell your operating system where to look for programs, so that you can just write `program` instead of `some_folders\\program`. But different operating systems have different ways to add a new directory to it. \n\nAdd the installation directory to PATH. To do so you must open again a terminal window as in Step 1. Run \n\n   ```\n   touch ~/.bash_profile; open ~/.bash_profile\n   ```\n\n   This will open a file known as `.bash_profile`, where you can copy and paste the following line  \n\n   ```\n   export PATH=\"$PATH:desired_installation_path/vchdfe/bin\"\n   ```\n\n   where `desired_installation_path` is the folder where you installed the executable in the previous steps. Next we save the file we just modified and close it. Finally, to make sure that it will load those changes without rebooting the computer, run the following line in the terminal\n\n   ```\n   source ~/.bash_profile\n   ```\n\n\n\n# Syntax of the program\n\nBefore this, make sure you've completed the installation steps shown above. To use the executable you only need to open the terminal. Windows users can press the keys Windows + R, then type \"powershell\" and press Enter. MacOS users can press COMMAND + SPACE to open spotlight search, and type terminal and hit RETURN. \n\nThe basic syntax of this command is \n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"vchdfe pathtodata [–option option_value]","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"where `[]` denote optional arguments. \n\nYou can see a list of examples provided in a section below or type in the powershell\n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"vchdfe –help","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"\nto see a complete list of available arguments:\n  \n    positional arguments:\n      path                  path to CSV file containing data\n\n    optional arguments:\n      --first_id FIRST_ID   column index in CSV file for the first ID\n                            The user should specify the most \"granular\"\n                            id here. For instance, in the AKM context,\n                            this id would correspond to worker ids.\n                            (type: Int64, default: 1)\n      --second_id SECOND_ID\n                            column index in CSV file for the second ID\n                            Use the less granular type (e.g. Firm in the \n                            AKM example)\n                            (type: Int64, default: 2)\n      --outcome_id OUTCOME_ID\n                            column index in CSV file for outcome (e.g.\n                            Wage). (type: Int64, default: 4)\n      --no_first_id_effects No computing and showing of first_id effects\n      --no_cov_effects      No computing and showing of covariance effects\n      --algorithm ALGORITHM\n                            type of algorithm: Exact or JLA. It defaults\n                            to be Exact if the number of observations is\n                            less than 5000, and JLA otherwise. (default:\n                            \"Default\")\n      --simulations SIMULATIONS\n                            number of simulations in the JLA algorithm. It\n                            defaults to 100 * log(#total fixed effect)\n                            (type: Int64, default: 0)\n      --header              CSV file contains header\n      --first_id_display FIRST_ID_DISPLAY\n                            The display text associated with first_id\n                            (e.g. Person). (default: \"Person\")\n      --second_id_display SECOND_ID_DISPLAY\n                            The display text associated with second_id\n                            (e.g. Firm) (default: \"Firm\")\n      --outcome_id_display OUTCOME_ID_DISPLAY\n                            The display text associated with outcome_id\n                            (e.g. Wage) (default: \"Wage\")\n      --detailed_output_path DETAILED_OUTPUT_PATH\n                            path to the CSV for the detailed output for\n                            each observable (default:\n                            \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\variance_components.csv\")\n      --results_path RESULTS_PATH\n                            path to the results of the output (default:\n                            \"C:\\\\Users\\\\owner\\\\Desktop\\\\vchdfe\\\\results.txt\")\n      --write_detailed_CSV  write the detailed output to a CSV\n      --write_results       write the results to a file\n      --print_level PRINT_LEVEL\n                            Level of verbosity of output. (type: Int64,\n                            default: 1)\n      -h, --help            show this help message and exit\n      \n      \n\nThe executable can provide two different type of outputs. The first (that we refer to as Results) is a log file, in txt format, that writes a summary of the main results : bias-corrected components and some statistics from the leave-out sample. The second type of output (that we refer to as Detailed Output) is a DataFrame, stored in csv format, that saves important objects computed throughout the algorithm such a subset of the original data corresponding to the leave-out sample, the corresponding Pii and Bii, and the fixed effects that can be used to compute the plug-in variance components. \n\nThe program provides two outputs: a file contains the main results and a detailed output file. \n\n## Results output\n\nBy using the argument `--write_results` followed by `--results_path MYFILE`, you can specify that you want to write the table with the main results in `MYFILE`. \n\n\n## Detailed output\n\nYou can save all the output details in a CSV file. You need to use `--write_detailed_CSV` followed by `--detailed_output_path MYFILE`, where `MYFILE` is the path to the detailed output file. \n\nThe detailed output file includes:                                                                                                                                                                                                                              \n                           \n- `observation` : observation identifier in the original dataset that belong to the Leave-out connected set. For instance, if the number 3 appears in this vector, it means that the third observation in the original dataset belongs to the leave-out connected set. \n- `first_id`: the first identifier corresponding to each observation in `obs` (e.g. worked ids in the leave-out connected set).\n- `second_id`: the second identifier corresponding to each observation in `obs` (e.g. firm ids in the leave-out connected set).\n- `D_alpha`: the fixed effect for the first identifier corresponding to each observation. \n- `F_psi`: the fixed effect for the second identifier corresponding to each observation. \n- `Pii`: statistical leverage corresponding to each observation in `obs`.\n- `Bii_first`: The Bii for the variance of `first_id` effects corresponding to each observation in `obs`.\n- `Bii_second`: The Bii for the variance of `second_id` effects corresponding to each observation in `obs`.\n- `Bii_cov`: The Bii for the co-variance of `first_id` and `second_id` effects corresponding to each observation in `obs`.\n\n# Typical Executable Workflow (Windows)\n\nYou begin by opening the Powershell (as administrator), and typing \n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"cd \"pathtodataset\"","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"\nThen, if you completed all instalation steps you may run\n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"vchdfe my_data.csv –OPTIONS","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"\nwhere `OPTIONS` depend on the structure of your data. \n\nIf you didn't complete step 5 of the installation, you will have to run instead \n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"installationfolder\\vchdfe\\bin\\vchdfe mydata.csv –OPTIONS","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"\n## Detailed examples\n\nSuppose we have a dataset `my_data.csv` that is stored in `\"project_path\"`. The structure of the data is such that the fifth column corresponds to the worker identifier and the third column corresponds to the firm identifier. Moreover, we want to store a summary of the results in this location `\"project_path/summary.txt\"`, and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  `\"project_path/output.csv\"`. \n\n\n1. To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the terminal \n\n   ```\n   cd project_path\n   \n   vchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt\n   ```\n\n2. To run the same thing while specifying 1000 simulations for the JLA algorithm that estimates (Pii,Bii) described in the computational appendix of KSS, we type in the terminal \n\n   ```\n   cd project_path\n\n   vchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt --algorithm JLA --simulations 1000\n   ```\n\n3. To only obtain the bias-correction for the variance of firm effects, we type in the powershell \n\n   ```\n   cd project_path\n\n   vchdfe my_data.csv --first_id 5 --second_id 3 --no_first_effects --no_cov_effects --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt  --algorithm JLA --simulations 1000\n   ```\n\n\n# Typical Executable Workflow (MacOS)\n\nYou begin by opening the Command Line (COMMAND + SPACE), and typing \n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"cd pathtodataset","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"\nThen, if you completed all instalation steps you may run\n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"vchdfe my_data.csv –OPTIONS","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"\nwhere `OPTIONS` depend on the structure of your data. \n\nIf you didn't complete step 4 of the installation, you will have to run instead \n","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"installationfolder/vchdfe/bin/vchdfe mydata.csv –OPTIONS ```","category":"page"},{"location":"Executable/#Detailed-examples","page":"Executable","title":"Detailed examples","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"Suppose we have a dataset my_data.csv that is stored in \"project_path\". The structure of the data is such that the fifth column corresponds to the worker identifier and the third column corresponds to the firm identifier. Moreover, we want to store a summary of the results in this location \"project_path/summary.txt\", and the detailed output (Pii and Bii matrices, stored fixed effects, etc.) here  \"project_path/output.csv\". ","category":"page"},{"location":"Executable/","page":"Executable","title":"Executable","text":"To obtain all three bias-corrected components (variance of worker effects, variance of firm effects, and covariance of worker-firm effects), we only need to type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt\nTo run the same thing while specifying 1000 simulations for the JLA algorithm that estimates (Pii,Bii) described in the computational appendix of KSS, we type in the terminal \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt --algorithm JLA --simulations 1000\nTo only obtain the bias-correction for the variance of firm effects, we type in the powershell \ncd project_path\n\nvchdfe my_data.csv --first_id 5 --second_id 3 --no_first_effects --no_cov_effects --write_results --write_detailed_CSV --detailed_output_path output.csv --results_path summary.txt  --algorithm JLA --simulations 1000","category":"page"},{"location":"Executable/#About-the-current-version","page":"Executable","title":"About the current version","text":"","category":"section"},{"location":"Executable/","page":"Executable","title":"Executable","text":"The bias-correction currently only runs on a two-way fixed effects model without controls.\nThe user can manually preadjust the outcome. For instance, in an AKM context, the user can run first   y = pe + fe + Xb + e , where pe are person effects and fe are firm effects, and feed into the routine y-Xb as the outcome.","category":"page"},{"location":"#VarianceComponentsHDFE","page":"Home","title":"VarianceComponentsHDFE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package estimates a two-way fixed effects model and computes its associated variance components using the methodology developed by Kline, Saggio and Sølvsten (KSS). We provide the usual Julia Package as well as an executable/app that can be run in the terminal. For more details about this please see the corresponding Executable section. The link to the repository can be found  here.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please report any bug or problem here.","category":"page"},{"location":"#Contents","page":"Home","title":"Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"Executable.md\", \"Package.md\"]\nDepth = 3","category":"page"},{"location":"Package/#Contents","page":"Package","title":"Contents","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"Pages = [\"Package.md\"]\nDepth = 3","category":"page"},{"location":"Package/#Setting-up-Development-Enviroment","page":"Package","title":"Setting up Development Enviroment","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"Install Julia and GitHub Desktop - not strictly required but never hurts to have it!\nInstall vscode and follow basic instructions in https://github.com/ubcecon/tutorials/blob/master/vscode.md\nIn particular, https://github.com/ubcecon/tutorials/blob/master/vscode.md#julia, making sure to do the code formatter step.\nand the git settings in https://github.com/ubcecon/tutorials/blob/master/vscode.md#general-packages-and-setup\nClone the repo by either:\nClicking on the Code then Open in GitHub Desktop.\nAlternatively, you can go ] dev https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl in a Julia REPL and it will clone it to the .julia/dev folder.\nIf you wanted to, you could then drag that folder back into github desktop.\nOpen it in vscode by right-clicking on the folder it installed to, and then opening a vscode project.\nOpen the Julia repl in vscode  (Ctrl-Shift-P and then go Julia REPL or something to find it.\ntype ] instantiate to install all of the packages.  Get coffee.\nIn the REPL run ] test and it should do the full unit test.","category":"page"},{"location":"Package/#Functions-in-this-package","page":"Package","title":"Functions in this package","text":"","category":"section"},{"location":"Package/#Main-Functions","page":"Package","title":"Main Functions","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"    leave_out_estimation(y,first_id,second_id,controls,settings)\n    get_leave_one_out_set(y, first_id, second_id, settings, controls)    ","category":"page"},{"location":"Package/#VarianceComponentsHDFE.leave_out_estimation-NTuple{5,Any}","page":"Package","title":"VarianceComponentsHDFE.leave_out_estimation","text":"leave_out_estimation(y, first_id, second_id, controls, settings)\n\n\nReturns the bias-corrected components, the vector of coefficients, the corresponding fixed effects for every observation, and the diagonal matrices containing the Pii and Biis. \n\nArguments\n\ny: outcome vector\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on VCHDFESettings\ncontrols: at this version only controls=nothing is supported.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.get_leave_one_out_set-NTuple{5,Any}","page":"Package","title":"VarianceComponentsHDFE.get_leave_one_out_set","text":"get_leave_one_out_set(y, first_id, second_id, settings, controls)\n\n\nReturns a tuple with the observation number of the original dataset that belongs to the Leave-out connected set as described in Kline,Saggio, Solvesten. It also provides the corresponding outcome and identifiers in this connected set. \n\nArguments\n\ny: outcome vector\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on VCHDFESettings\ncontrols: at this version only controls=nothing is supported.\n\n\n\n\n\n","category":"method"},{"location":"Package/#Auxiliary-Functions","page":"Package","title":"Auxiliary Functions","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"    find_connected_set(y, first_idvar, second_idvar, settings)\n    prunning_connected_set(yvec, first_idvar, second_idvar, obs_id, settings)\n    drop_single_obs(yvec, first_idvar, second_idvar,obs_id)\n    compute_movers(first_id,second_id)\n    eff_res(::ExactAlgorithm, X,first_id,second_id,match_id, K, settings)\n    eff_res(lev::JLAAlgorithm, X,first_id,second_id,match_id, K, settings)\n    compute_matchid(second_id,first_id)   ","category":"page"},{"location":"Package/#VarianceComponentsHDFE.find_connected_set-NTuple{4,Any}","page":"Package","title":"VarianceComponentsHDFE.find_connected_set","text":"find_connected_set(y, first_idvar, second_idvar, settings)\n\n\nReturns a tuple of observation belonging to the largest connected set with the corresponding identifiers and outcomes. This requires to have the data sorted by first identifier, and time period (e.g. we sort by worked id and year). This is also the set where we can run AKM models with the original data.\n\nArguments\n\ny: outcome (e.g. log wage)\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.prunning_connected_set-NTuple{5,Any}","page":"Package","title":"VarianceComponentsHDFE.prunning_connected_set","text":"prunning_connected_set(yvec, first_idvar, second_idvar, obs_id, settings)\n\n\nThis function prunes the dataset from articulation points. If the first identifier is worker id it means that it prunes workers that would disconnect the graph if they were dropped.\n\nArguments\n\nyvec: outcome (e.g. log wage)\nfirst_idvar: first identifier (e.g. worker id)\nsecond_idvar: second identifier (e.g. firm id)\nobs_id: observation identifier.\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.drop_single_obs-NTuple{4,Any}","page":"Package","title":"VarianceComponentsHDFE.drop_single_obs","text":"drop_single_obs(yvec, first_idvar, second_idvar, obs_id)\n\n\nThis function drops observations that correspond with first identifiers with a single observation. For example, if first identifier is worker id, it will drop observations for workers that only appear once in the data.\n\nArguments\n\nyvec: outcome (e.g. log wage)\nfirst_idvar: first identifier (e.g. worker id)\nsecond_idvar: second identifier (e.g. firm id)\nobs_id: observation identifier.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.compute_movers-Tuple{Any,Any}","page":"Package","title":"VarianceComponentsHDFE.compute_movers","text":"compute_movers(first_id, second_id)\n\n\nReturns a vector that indicates whether the first_id (e.g. worker) is a mover across second_id (e.g. firms), as well as a vector with the number of periods that each first_id appears.\n\nArguments\n\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.eff_res-Tuple{ExactAlgorithm,Any,Any,Any,Any,Any,Any}","page":"Package","title":"VarianceComponentsHDFE.eff_res","text":"eff_res(_, X, first_id, second_id, match_id, K, settings)\n\n\nThis function computes the diagonal matrices containing Pii and Bii under the Exact Algorithm. See appendix in KSS for more information.\n\nArguments\n\nX: the design matrix in the linear model.\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nmatch_id: match identifier between first and second identifier.For example, match can be the identifier of every worker-firm combination.\nK: number of covariates in addition to the fixed effects. Currently only 0 is supported.\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.eff_res-Tuple{JLAAlgorithm,Any,Any,Any,Any,Any,Any}","page":"Package","title":"VarianceComponentsHDFE.eff_res","text":"eff_res(lev, X, first_id, second_id, match_id, K, settings)\n\n\nThis function computes the diagonal matrices containing Pii and Bii under Johnson-Linderstrauss Algorithm. See appendix in KSS for more information.\n\nArguments\n\nX: the design matrix in the linear model.\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\nmatch_id: match identifier between first and second identifier.For example, match can be the identifier of every worker-firm combination.\nK: number of covariates in addition to the fixed effects. Currently only 0 is supported.\nsettings: settings based on data type VCHDFESettings. Please see the reference provided below.\n\n\n\n\n\n","category":"method"},{"location":"Package/#VarianceComponentsHDFE.compute_matchid-Tuple{Any,Any}","page":"Package","title":"VarianceComponentsHDFE.compute_matchid","text":"compute_matchid(second_id, first_id)\n\n\nComputes a match identifier for every combination of first and second identifier. For example, this can be the match identifier of worker-firm combinations.\n\nArguments\n\nfirst_id: first identifier (e.g. worker id)\nsecond_id: second identifier (e.g. firm id)\n\n\n\n\n\n","category":"method"},{"location":"Package/#Datatypes-in-this-package","page":"Package","title":"Datatypes in this package","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"    ExactAlgorithm\n    JLAAlgorithm\n    VCHDFESettings\n    ","category":"page"},{"location":"Package/#VarianceComponentsHDFE.ExactAlgorithm","page":"Package","title":"VarianceComponentsHDFE.ExactAlgorithm","text":"struct ExactAlgorithm <: AbstractLeverageAlgorithm\n\nData type to pass to VCHDFESettings type, to indicate Exact algorithm\n\n\n\n\n\n","category":"type"},{"location":"Package/#VarianceComponentsHDFE.JLAAlgorithm","page":"Package","title":"VarianceComponentsHDFE.JLAAlgorithm","text":"struct JLAAlgorithm <: AbstractLeverageAlgorithm\n\nData type to pass to VCHDFESettings type, to indicate JLA algorithm\n\nFields\n\nnum_simulations: number of simulations in estimation. If num_simulations = 0, defaults to 100 * log(#total fixed effect)\"\n\n\n\n\n\n","category":"type"},{"location":"Package/#VarianceComponentsHDFE.VCHDFESettings","page":"Package","title":"VarianceComponentsHDFE.VCHDFESettings","text":"struct VCHDFESettings{LeverageAlgorithm}\n\nThe VCHDFESettings type is to pass information to methods regarding which algorithm to use. \n\nFields\n\ncg_maxiter: maximum number of iterations (default = 300)\nleverage_algorithm: which type of algorithm to use (default = ExactAlgorithm())\nfirst_id_effects: includes first id effects. At this version it is required to include the firstideffects. (default = true)\ncov_effects: includes covariance of first-second id effects. At this version it is required to include the cov_effects. (default = true)\nprint_level: prints the state of the program in std output. If print_level = 0, the app prints nothing in the std output. (default = 1)\nfirst_id_display_small: name of the first id in lower cases (default = person)\nfirst_id_display: name of the first id (default = Person)\nsecond_id_display_small: name of the second id in lower cases (default = firm)\nsecond_id_display: name of the second id (default = Firm)\noutcome_id_display_small: name of the observation id in lower cases (default = wage)\noutcome_id_display: name of the observation id (default = Wage)\n\n\n\n\n\n","category":"type"},{"location":"Package/#Typical-Julia-Workflow","page":"Package","title":"Typical Julia Workflow","text":"","category":"section"},{"location":"Package/","page":"Package","title":"Package","text":"#Load the required packages\nusing VarianceComponentsHDFE, DataFrames, CSV\n\n#Load dataset\ndata = DataFrame!(CSV.File(\"dataset.csv\"; header=false))\n\n#Extract vectors of outcome, workerid, firmid\nid = data[:,1]\nfirmid = data[:,2]\ny = data[:,3]\n\n#Define the settings using our structure: JLA Algorithm\nsettings_JLA = VCHDFESettings(leverage_algorithm = JLAAlgorithm(), first_id_effects=true, cov_effects=true)\n\n#Define the settings using our structure: Exact Algorithm\nsettings_Exact = VCHDFESettings(leverage_algorithm = ExactAlgorithm(), first_id_effects=true, cov_effects=true)\n\n#Compute Leave-Out Connected Set\nobs,  y , id , firmid, controls = get_leave_one_out_set(y, id, firmid, settings_JLA, nothing)\n\n#Run Leave-Out Correction in the Leave-Out Set\nθ_worker, θ_firms, θCOV, β, Dalpha, Fpsi, Pii, Bii_first, Bii_second, Bii_cov = leave_out_estimation(y,id,firmid,nothing,settings_JLA)\n\n#Print the Bias-Corrected Components from the output\nprintln(\"Bias-corrected Variance of Worker Effects:\", θ_worker)\nprintln(\"Bias-corrected Variance of Firm Effects:\", θ_firms)\nprintln(\"Bias-corrected Variance of Cov Worker-Firm Effects:\", θCOV)\n","category":"page"}]
}

# Development Environment Setup

## Setting up Development Enviroment
1. Install [Julia](https://julialang.org/downloads/) and [GitHub Desktop](https://desktop.github.com/) - not strictly required but never hurts to have it!
2. Install `vscode` and follow basic instructions in https://github.com/ubcecon/tutorials/blob/master/vscode.md
  - In particular, https://github.com/ubcecon/tutorials/blob/master/vscode.md#julia, making sure to do the code formatter step.
  - and the git settings in https://github.com/ubcecon/tutorials/blob/master/vscode.md#general-packages-and-setup
3. Clone the repo by either:
  - Clicking on the `Code` then `Open in GitHub Desktop`.
  - Alternatively, you can go `] dev https://github.com/HighDimensionalEconLab/VarianceComponentsHDFE.jl` in a Julia REPL and it will clone it to the `.julia/dev` folder.
  - If you wanted to, you could then drag that folder back into github desktop.
4. Open it in vscode by right-clicking on the folder it installed to, and then opening a vscode project.
5. Open the Julia repl in vscode  (`Ctrl-Shift-P` and then go `Julia REPL` or something to find it.
6. type `] instantiate` to install all of the packages.  Get coffee.
6. In the REPL run `] test` and it should do the full unit test.

## Formatting code
- Assuming that you setup the code formatter correctly in the above instructions, before checking in any code you should go `Ctrl-Shift-P` and type `Formatter`.

## Code Standards and Design Principles
- Use the unicode math symbol matching the algebra whenever possible.
- Follow https://github.com/jrevels/YASGuide largely.
    - https://github.com/QuantEcon/lecture-source-jl/blob/master/style.md also is useful, but is intended more for "scripting" code rather than package code.


## Code Style
- Use the unicode and/or math symbol matching the algebra whenever possible
- Follow https://github.com/jrevels/YASGuide where possible
    - https://github.com/QuantEcon/lecture-source-jl/blob/master/style.md also is useful, but is intended more for "scripting" code rather than package code.
- The `.JuliaFormatter.toml` file has settings for the automatic formatting of the code
  - **CURRENTLY BROKEN**
  - To use it in Atom/Juno, use `<Ctrl+P>` then type in `format` and you will find it.
  - Before formatting code, consider ensureing the unit tests pass.  I haven't seen bugs with the formatting yet, but there may be some.
  - For vscode, install https://marketplace.visualstudio.com/items?itemName=singularitti.vscode-julia-formatter then in your    `settings.json` add in (following the instructions there)
    ```
      "[julia]": {
        "editor.defaultFormatter": "singularitti.vscode-julia-formatter"
    },
    ```
  - In vscode if you then go `ctrl-p` and type `format` you will be able to choose the singularitti one as your default julia editor
- It is essential to force the use of `LF` for anyone using Windows on any computer, to do this
  - Git config in https://julia.quantecon.org/more_julia/version_control.html#Setup
  - Atom settings in https://julia.quantecon.org/more_julia/tools_editors.html#Installing-Atom
  - VSCode settings (i.e. files.eof) in https://github.com/ubcecon/tutorials/blob/master/vscode.md#general-packages-and-setup

  ## Compiling an Executable
  - To compile an executable of the package run the following Julia commands.
  ```julia
    using PackageCompiler
    create_app("PathToPackageDir/VarianceComponentsHDFE", "PathToExecutableDir/VarianceComponentsHDFEExecutable)
  ```
  

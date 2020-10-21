Notes:

* First benchmark is slow because it needs to create a `tune.json` file.

* The preamble for all of the following code is:

```
using PkgBenchmark, VarianceComponentsHDFE
```

* To do it for a specific commit or version, do something like the following (**note** the benchmarks must have existed at that commit, as it will do a checkout):

```
benchmarkpkg(VarianceComponentsHDFE, "commit-hash")
```

* To compare two commits:

```
judge(VarianceComponentsHDFE, "new-commit-hash", "old-commit-hash")
```

First argument is optional; will use state of repo if only one hash given.

Can also provide an `f` like median, etc. as a kwarg, for the judging function.

Gives output as a `Benchmarkjudement` which can also be written to Markdown:

* To get Markdown, use something like this:

```
export_markdown("benchmarks/trial.md", benchmarkpkg(VarianceComponentsHDFE))
```
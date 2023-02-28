# Architecture for Linken

The previous pipelines  composed of `Makefiles` and shell scripts for preparing and running the analyses, followed by Perl scripts to generate data in format that can be used for plotting and statistical analyses.
To be clear, the previous pipelines had no issues.
Linken was written to improve the existing pipeline in the following aspects: (1) encapsulation through containerization, (2) cleaner API, and (3) improving speed.

The first aspect, encapsulation through containerization, is critical to ensure safety and robustness of the pipeline environment.
The previous pipeline used shell environment variables (`env`), which is ephemeral.
Through containerization, the variables are instead statically linked, therefore reducing the likelihood of a bug that could occur through improper `env` declaration.

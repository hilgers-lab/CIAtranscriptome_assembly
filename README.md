# CIA Transcriptome Assembly

Snakefile pipeline of all the steps taken to reconstruct the CIA transcriptome assembly as in (Alfonso-Gonzalez, 2022). 

# Documentation

The documentation for this project is hosted at [ReadTheDocs](https://cia-transcriptome-assembly.readthedocs.io/en/latest/).

# Troubleshooting

Occasionally there can be an issue installing R packages in the `cia-sqanti`
environment. This will manifest in an error like this: 

```
Error in if (nzchar(SHLIB_LIBADD)) SHLIB_LIBADD else character() : 
argument is of length zero
```

if you run into this error, follow the instructions from [this thread](https://stackoverflow.com/questions/53813323/installing-r-packages-in-macos-mojave-error-in-if-nzcharshlib-libadd/54778735?stw=2#54778735).

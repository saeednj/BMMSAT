# MaxSAT Evaluation 2018 Patch

This patch was done by Ruben Martins to:
- catch SIGTERM signal 
- do not print intermediate solutions until SIGTERM is received
- faster printing of model using stringstreams

To patch SATLike execute the following command:

```
patch -p1 < SATLike-SIGTERM.patch
```

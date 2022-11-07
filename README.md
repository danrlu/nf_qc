# qc
Quick QC and trim for raw sequencing reads

```
nextflow run danrlu/qc --se --input_folder /full_path/data
```

default is paired-end, so only single end needs `--se`

output will be in ${params.input_folder}/trimmed

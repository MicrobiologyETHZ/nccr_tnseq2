import pandas as pd




rule filter_data:
    input: cnts = OUTDIR/'counts/{sample}_counts.csv',
         cntrl =  OUTDIR/'counts/{sample}_counts.control.csv'
    output:
    params: ''
    shell: "python analysis.py filter {input.cnts} {input.cntrl} {output.sampleData} {output.exprData} "


rule runDEseq:
    input: sampleData = 'cnt_file',
         exprData = 'cntrl_file'
    output: sampleData = temp(sdf),
          exprData = temp(edf)
    params: ''
    shell: "python analysis.py compare {input.cnts} {input.cntrl} {output.sampleData} {output.exprData} "


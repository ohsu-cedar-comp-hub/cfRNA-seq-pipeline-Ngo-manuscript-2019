import pandas as pd

files = snakemake.input
out = snakemake.output[0]

def split_exon_fraction(files, out):
    sample_list = []
    for file in files:
        sample = file.split('.')[0]
        name = sample.split('/')[3]
        results = {}
        equal = False
        for line in open(file):
            if '==' in line:
                equal = True

            if not equal:
                line = line.split()
                key ='_'.join(line[:2])
                results[key] = float(line[-1])

            if equal:
                if not '==' in line and not'Group' in line:
                    line = line.split()
                    key = line[0]
                    results[key] = float(line[2])
                else:
                    continue
        Total = results['CDS_Exons'] + results["5'UTR_Exons"] + results["3'UTR_Exons"] + results['Introns'] + results['TSS_up_1kb'] + results['TSS_up_5kb'] + results['TSS_up_10kb'] + results['TES_down_1kb'] + results['TES_down_5kb'] + results['TES_down_10kb']
        exon_fraction = (results['CDS_Exons'] + results["5'UTR_Exons"] + results["3'UTR_Exons"]) / Total
        intron_fraction = results['Introns'] / Total
        intergenic_fraction = (results['TSS_up_1kb'] + results['TSS_up_5kb'] + results['TSS_up_10kb'] + results['TES_down_1kb'] + results['TES_down_5kb'] + results['TES_down_10kb']) / Total
        sample_list.append([name, exon_fraction, intron_fraction, intergenic_fraction])

    frame = pd.DataFrame(sample_list)
    frame.columns = ['Sample','Exon','Intron','Intergenic']
    frame.to_csv(out,sep='\t',index=False)

split_exon_fraction(files, out)

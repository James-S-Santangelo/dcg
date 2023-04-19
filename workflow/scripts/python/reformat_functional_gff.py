# Script to reformat GFF file

import re
import gffutils

# Creat dictionary mapping products to EC numbers
EC_num_to_products = {}
with open(snakemake.input['ec'][0], 'r') as fin:
    lines = fin.readlines()
    for l in lines:
        if l.startswith('ID'):
            EC_number = l.split('   ')[1].strip()
        if l.startswith('DE'):
            prod = l.split('   ')[1].strip()
            EC_num_to_products[EC_number] = prod
        else:
            pass

locus_tag_prefix = snakemake.params['locus_tag']
with open(snakemake.output[0], 'w') as fout:
    fout.write("##gff-version 3\n")
    gene_num = 1
    for feat in gffutils.DataIterator(snakemake.input['gff']):
        if feat.featuretype == 'gene':
            # Set CDS counter
            cds_count = 1
            
            # Add locus tag to genes
            feat.attributes["locus_tag"] = f'{locus_tag_prefix}_' + str(int(gene_num)).zfill(5)
            gene_num += 1
        elif feat.featuretype == 'mRNA':
            if 'product' in feat.attributes.keys() and 'EC_number' in feat.attributes.keys():
                # Reassign product if EC Number to 4 digits, else delete EC Number and keep hypothetical
                EC_number = feat.attributes['EC_number'][0]
                prod = feat.attributes['product'][0]
                if prod == 'hypothetical protein':
                    EC_split = EC_number.split('.')
                    if len(EC_split) < 4:
                        del feat.attributes['EC_number']
                    elif len(EC_split) == 4:
                        new_product = EC_num_to_products[EC_number]
                        # Handle cases where EC numbers have been reassigned
                        if 'Transferred' in new_product:
                            pattern = r'(?<=:\s)(\d\.\d+\.\d+\.\d+(?=\.))'
                            new_ec = re.search(pattern, new_product).group(1)
                            new_product = EC_num_to_products[new_ec]
                        new_product = new_product.replace('.', '')
                        feat.attributes['product'] = [new_product]
                        
                        # Add note to GFF about origin of Product
                        feat.attributes['note'] = ["Funnanotate product changed from hypothetical protein based on EC_number"]
        elif feat.featuretype == 'CDS':
            # Assign increment CDS to ensure unique ID
            feat.attributes['ID'][0] = re.sub('cds', f'cds{cds_count}', feat.attributes['ID'][0])
            cds_count += 1
        else:
            pass
        fout.write(str(feat) + '\n')

# Python script to remap attributes and fix some issues flagged by NCBI

import re
import gffutils

db = gffutils.create_db(snakemake.input['gff'][0], snakemake.output['db'], force = True)


gene_name_dict = {
    'PROFILIN4' : 'PFN4',
    'PEPTIDEN4(NACETYLBETAGLUCOSAMINYL)ASPARAGINE' : 'PNG1',
    'SEDOHEPTULOSE1' : 'SBPASE',
    'DELTA1PYRROLINE5CARBOXYLATE' : 'P5CS1',
    'PROFILIN1' : 'PFN1',
    'PROFILIN2' : 'PFN2',
    '3ISOPROPYLMALATE' : 'IMDH'
}

def fix_product(feature):
    
    # Get product based on feature type
    featuretype = feature.featuretype
    # mRNA has product as attribute
    if featuretype == 'mRNA':
        product = ','.join(feat.attributes['product'])
    # Have to get product from CDS's mRNA parent
    elif featuretype == 'CDS':
        mrna_parent = [f.id for f in db.parents(feat.id, featuretype='mRNA')][0]
        product = ','.join(db[mrna_parent].attributes['product'])

    # Get gene parent, which is required for some specific fixes
    gene_parent = [f.id for f in db.parents(feat.id, featuretype='gene')][0]
        
    # Some simple string substitutions for common problems
    if re.search('\sGRAS$', product):
        new_product = 'Scarecrow-like 8'
    elif 'Belongs' in product:
        new_product = re.sub('Belongs to the ', '', product)
    elif 'proteins' in product:
        new_product = re.sub('(?<=\s)proteins$', 'protein', product)
    elif 'domains' in product:
        new_product = re.sub('(?<=\s)domains$', 'domain', product)
    elif 'AT3g18370/MYF24_8' in product:
        new_product = re.sub('AT3g18370/MYF24_8', 'C2 domain-containing protein', product)
    elif re.search('\sA[T|t][a-z0-9]{7}', product):
        new_product = re.sub('\sA[T|t][a-z0-9]{7}', '', product)
    elif 'YGL059W' in product:
        new_product = re.sub(' YGL059W', '', product)
    elif product == 'glutathione transferase activity':
        new_product = 'glutathione transferase'
    elif product == 'negative regulation of ubiquitin-protein transferase activity':
        new_product = 'uv-b-insensitive 4'
    elif product == 'Control of topological states of DNA by transient breakage and subsequent rejoining of DNA strands':
        new_product = 'DNA topoisomerase 2'

    elif 'and other' in product:
        new_product = re.sub('\(and other\)\s', '', product)
    elif re.search('Pseudo\s', product):
        new_product = re.sub('Pseudo\s', '', product)
    elif product == 'Short calmodulin-binding motif containing conserved Ile and Gln residues.':
        new_product = 'IQ-domain 7'
    elif product == 'Viral Genome polyprotein':
        new_product = 'asparagine--tRNA ligase'
    elif 'HOMOLOG' in product:
        new_product = re.sub('\sHOMOLOG', '', product)
    elif 'Allatostatins [Cleaved' in product:
        new_product = 'anthranilate synthase 2'
    elif 'ditrans,polycis-polyprenyl diphosphate synthase' in product:
        new_product = 'ditrans,polycis-polyprenyl diphosphate synthase'
    
    # Problems with specific product annotations
    elif gene_parent == 'ACLI19_g3593':
        new_product = re.sub('\sactivity', '', product)
    elif gene_parent in ['ACLI19_g1358', 'ACLI19_g44008']:
        new_product = 'DNA topoisomerase 6 subunit A'
    elif gene_parent =='ACLI19_g20691':
        new_product = 'Putative SNAP25 homologous protein SNAP33'
    elif gene_parent == 'ACLI19_g43419':
        new_product = 'Putative SNAP25 homologous protein SNAP30'
    elif gene_parent == 'ACLI19_g47730':
        new_product = 'alpha-1 tubulin'
    else:
        new_product = product

    # Update product if present, create if not
    if 'product' in feature.attributes.keys():
        feat.attributes['product'][0] = new_product
    else:
        feat.attributes['product'] = [new_product]
        
    return(feat)

with open(snakemake.output['gff'], 'w') as fout:
    fout.write("##gff-version 3\n")
    for feat in db.all_features():
        
        if feat.featuretype == 'gene':
            # Remove incremental numbering of gene names generated by funannotate
            # Reassign as "gene" attribute
            if 'Name' in feat.attributes.keys():
                feat.attributes['gene'] = [re.sub('_\d+$', '', feat.attributes['Name'][0])]
                feat.attributes['gene'][0] = feat.attributes['gene'][0].upper()
                feat.attributes.pop('Name')
                
                # Fix verbose/incorect gene names
                gene_name_base = feat.attributes['gene'][0].split('_')[0]
                if gene_name_base in gene_name_dict.keys():
                    feat.attributes['gene'][0] = re.sub(re.escape(gene_name_base), 
                                                        gene_name_dict[gene_name_base], 
                                                        feat.attributes['gene'][0])

                
            # Extend truncated gene to stop codon
            if feat.id == 'ACLI19_g45676':
                feat.start = 1814
            # Mark terminal gene with missing stop codong as "pseudo"
            if feat.id == 'ACLI19_g5041q1':
                feat.attributes['pseudo'] = ['true']            
                
        elif feat.featuretype == 'mRNA':
            # Remove first transcript of ACLI19_g11316 since missing stop codon and second isoform looks good
            if feat.id == 'ACLI19_g11316.t1':
                continue
                
            # Assign transcript IDs
            locus_tag = [f.attributes['locus_tag'][0] for f in db.parents(feat.id, featuretype='gene')][0]
            transcript = feat.attributes['ID'][0].split('.')[1]
            feat.attributes['transcript_id'] = [f"gnl|NessUTM|mrna{transcript}.{locus_tag}"]

            # Extend truncated mRNA to stop codon
            if feat.id == 'ACLI19_g45676.t1':
                feat.start = 1814
            
            # Fix and reassign product
            feat = fix_product(feat)
            
        elif feat.featuretype == 'CDS':
            
            # Remove first transcript of ACLI19_g11316 since missing stop codon and second isoform looks good
            if feat.attributes['Parent'][0] == 'ACLI19_g11316.t1':
                continue
                
            # Extend truncated CDS to stop codon
            if feat.id == 'ACLI19_g45676.t1.cds2':
                feat.start = 1814
                
            # Mark terminal gene missing stop codong as "pseudo"
            gene_parent = [f.id for f in db.parents(feat.id, featuretype='gene')][0]
            if gene_parent == 'ACLI19_g50411':
                feat.attributes['pseudo'] = ['true']
                
            # Assign protein IDs
            locus_tag = [f.attributes['locus_tag'][0] for f in db.parents(feat.id, featuretype='gene')][0]
            transcript = feat.attributes['Parent'][0].split('.')[1]
            feat.attributes['protein_id'] = [f"gnl|NessUTM|{transcript}.{locus_tag}"]
            
            # Get, fix, and assign product
            feat = fix_product(feat)
            
            # Assign EC Number to CDS features
            mrna_parent = [f.id for f in db.parents(feat.id, featuretype='mRNA')][0]
            if 'ec_number' in db[mrna_parent].attributes.keys():
                feat.attributes['ec_number'] = db[mrna_parent].attributes['ec_number']

        # Don't need exons in final GFF
        elif feat.featuretype == 'exon':
            continue
        
        fout.write(f"{str(feat)}\n")

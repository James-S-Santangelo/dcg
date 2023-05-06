import re
import gffutils

db = gffutils.create_db(snakemake.input['gff'], snakemake.output['db'], force = True)

genes_to_remove = []
# Iterate over gene features
for gene in db.features_of_type('gene'):
    
    # Get all genes completely overlapping this gene's coordinates
    # Note: This will return the gene itself
    overlapping_genes = list(db.features_of_type(strand = gene.strand, 
                                                 limit = f"{gene.seqid}:{gene.start}-{gene.end}", 
                                                 featuretype='gene', 
                                                 completely_within=True))
    
    # If there is more than 1 gene completely overlapping this span
    if len(overlapping_genes) > 1:
        # Iterate over the genes overlapping this span
        over_gene_length_dict = {}
        for over_gene in overlapping_genes:
            
            # Get length of all genes overlapping span
            gene_length = over_gene.end - over_gene.start
            over_gene_length_dict[over_gene.id] = gene_length
        
        # Let's only keeps the longest gene
        # This works because the nested genes are always smaller
        gene_toKeep = max(over_gene_length_dict, key=over_gene_length_dict.get)
        nested_genes_to_consider = [k for k in over_gene_length_dict.keys() if k != gene_toKeep]

        # Iterate again but only over small nested genes. Count CDS and check CDS overlap
        for over_gene in overlapping_genes:
            if over_gene.id in nested_genes_to_consider:
                
                # Get list with all CDSs of nested gene
                cds_list = [cds for cds in db.children(over_gene.id, featuretype='CDS')]
                
                # Iterate over Nested CDSs
                any_cds_overlap = 0
                for cds in cds_list:
                    
                    # Get CDS features overlapping nested CDS
                    overlapping_cds = list(db.features_of_type(strand = cds.strand, 
                                                               limit = f"{cds.seqid}:{cds.start}-{cds.end}", 
                                                               featuretype = 'CDS', 
                                                               completely_within=False))
                    
                    # If the CDS has other overlapping CDS, flag for removal
                    # Uses flag to ensure all CDSs considered (i.e., not all may overlap)
                    if len(overlapping_cds) >= 2:
                        any_cds_overlap = 1
                    else:
                        continue
                
                # Add gene to removal list if it overlaps existing CDS
                # Only add if not already there. Handles cases of multiple overlapping nested genes
                if any_cds_overlap == 1:
                    if over_gene.id not in genes_to_remove:
                        genes_to_remove.append(over_gene.id)
                else:
                    continue

with open(snakemake.output['gff'], 'w') as fout:
    fout.write("##gff-version 3\n")
    gene_num = 1
    for feat in db.all_features():
        
        if feat.featuretype == 'gene':
            
            # Remove gene features
            if feat.id in genes_to_remove:
                continue
            else:
                # Re-map ID, Locus Tag, and write
                new_id = f"ACLI19_g{gene_num}"
                new_locus_tag = f'P8452_' + str(int(gene_num)).zfill(5)
                feat.attributes['ID'][0] = new_id
                feat.attributes['locus_tag'][0] = new_locus_tag
                fout.write(f"{str(feat)}\n")
                gene_num += 1

        elif feat.featuretype == 'mRNA':
            
            # Remove features if gene parent is removed
            gene_parent = [f.id for f in db.parents(feat.id, featuretype='gene')][0]            
            if gene_parent in genes_to_remove:
                continue
            else:
                feat.attributes['ID'][0] = re.sub("ACLI19_g\d+", new_id, feat.attributes['ID'][0])
                feat.attributes['Parent'][0] = re.sub("ACLI19_g\d+", new_id, feat.attributes['Parent'][0])
                feat.attributes['transcript_id'][0] = re.sub("P8452_\d{5}", new_locus_tag, feat.attributes['transcript_id'][0])
                fout.write(f"{str(feat)}\n")
                
        elif feat.featuretype == 'CDS':
            
            # Remove features if gene parent is removed
            gene_parent = [f.id for f in db.parents(feat.id, featuretype='gene')][0]            
            if gene_parent in genes_to_remove:
                continue
            else:
                feat.attributes['ID'][0] = re.sub("ACLI19_g\d+", new_id, feat.attributes['ID'][0])
                feat.attributes['Parent'][0] = re.sub("ACLI19_g\d+", new_id, feat.attributes['Parent'][0])
                feat.attributes['protein_id'][0] = re.sub("P8452_\d{5}", new_locus_tag, feat.attributes['protein_id'][0])
                fout.write(f"{str(feat)}\n")



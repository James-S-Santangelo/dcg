import os
import gffutils

db = gffutils.create_db(snakemake.input['gff'][0],
                        './tmp.db',
                        force=True)

transcripts_to_toss = []
for gene in db.features_of_type('gene'):
    child_mrna = list(db.children(gene.id, featuretype='mRNA'))
    if len(child_mrna) > 1:
        isoform_lengths = {}
        for mrna in child_mrna:
            isoform_lengths[mrna.id] = sum([cds.end - cds.start for cds in db.children(mrna.id, featuretype='CDS')])
        if len(set(isoform_lengths.values())) == 1:
            same_length_to_toss = list(isoform_lengths.keys())[1:]
            transcripts_to_toss += same_length_to_toss
        elif len(set(isoform_lengths.values())) != 1:
            longest_isoform_length = max(isoform_lengths.values())
            to_toss = [k for k, v in isoform_lengths.items() if v < longest_isoform_length]
            transcripts_to_toss += to_toss
            id_match_max = [k for k, v in isoform_lengths.items() if v == longest_isoform_length]
            if len(id_match_max) > 1:
                transcripts_to_toss += id_match_max[1:]
    else:
        continue

with open(snakemake.output[0], 'w') as fout:
    fout.write("##gff-version 3\n")
    for feat in db.all_features():
        if feat.featuretype == 'mRNA':
            if feat.id in transcripts_to_toss:
                pass
            else:
                fout.write(f"{str(feat)}\n")
        elif feat.featuretype == 'CDS':
            mrna_parent = [f.id for f in db.parents(feat.id, featuretype='mRNA')][0]
            if mrna_parent in transcripts_to_toss:
                pass
            else:
                fout.write(f"{str(feat)}\n")
        else:
            fout.write(f"{str(feat)}\n")


os.remove('./tmp.db')

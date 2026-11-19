#!/bin/bash
O=/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119
UA='research-client/1.0'
g(){ mkdir -p "$(dirname "$2")"
  for i in 1 2 3 4 5 6; do curl -fsSL -C - --retry 3 --retry-delay 10 --connect-timeout 30 -A "$UA" -o "$2" "$1" && { echo "ok $2 $(stat -c %s "$2")"; return 0; }; sleep 20; done; echo "FAILED $1"; return 1; }
R=$O/raw
C=https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest
G=https://www.guidetopharmacology.org/DATA
P=https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data
g $C/chembl_37_sqlite.tar.gz $R/chembl/chembl_37_sqlite.tar.gz &
(for f in interactions.csv ligands.csv targets_and_families.csv ligand_id_mapping.csv; do g $G/$f $R/gtopdb/$f; done) &
(g https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv $R/complexportal/9606.tsv
 for f in interaction_input.csv complex_input.csv gene_input.csv protein_input.csv; do g $P/$f $R/cellphonedb/$f; done
 curl -fsSL -A "$UA" https://api.github.com/repos/ventolab/cellphonedb-data/commits/master -o $R/cellphonedb/master_commit.json
 g https://unmtid-dbs.net/download/DrugCentral/2027/drug.target.interaction.tsv.gz $R/drugcentral/drug.target.interaction.tsv.gz
 curl -fsSL -A "$UA" https://unmtid-dbs.net/download/DrugCentral/2027/ -o $R/drugcentral/directory_listing_2027.html
 g 'https://omnipathdb.org/interactions?genesymbols=1&datasets=ligrecextra,omnipath&fields=sources,references,datasets,type&license=academic' $R/omnipath/interactions.tsv
 g 'https://omnipathdb.org/intercell?scope=generic&aspect=functional&categories=ligand,receptor&license=academic' $R/omnipath/intercell_generic.tsv
 g 'https://omnipathdb.org/complexes?genesymbols=1&fields=sources,references&license=academic' $R/omnipath/complexes.tsv) &
wait
echo ALLDONE

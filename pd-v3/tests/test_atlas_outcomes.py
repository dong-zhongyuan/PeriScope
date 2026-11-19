import unittest
import numpy as np
import pandas as pd
from pdv3.atlas_outcomes import restore_original,normalize_ebi,normalize_native


class AtlasOrientationTests(unittest.TestCase):
    def native_row(self):
        return pd.DataFrame([dict(chromosome=23,base_pair_location=12345,
            effect_allele='T',other_allele='A',beta=-.31,standard_error=.04,
            p_value=1e-12,neg_log10_p_value=12.,variant_id='rs1',rsid='rs1',
            phenotype='p',source_y='Oxford',build_y='GRCh37',n_y=np.nan,
            n_y_status='X_not_assigned')])

    def test_native_schema_keeps_direction_coordinates_and_unknown_n(self):
        source=self.native_row();before=source.copy(deep=True)
        row=normalize_native(source).iloc[0]
        self.assertEqual((row.chrom_y,row.pos_y,row.build_y),('X',12345,'GRCh37'))
        self.assertEqual((row.ea_y,row.oa_y,row.beta_y,row.se_y),('T','A',-.31,.04))
        self.assertEqual((row.p_y,row.neg_log10_p_y,row.variant_y),(1e-12,12.,'rs1'))
        self.assertTrue(pd.isna(row.n_y));self.assertEqual(row.n_y_status,'X_not_assigned')
        pd.testing.assert_frame_equal(source,before)

    def test_native_rejects_wrong_build_and_ambiguous_columns(self):
        source=self.native_row();source['build_y']='GRCh38'
        with self.assertRaises(ValueError):normalize_native(source)
        source=self.native_row();source['beta_y']=.31
        with self.assertRaises(ValueError):normalize_native(source)
        with self.assertRaises(ValueError):normalize_native(self.native_row().drop(columns='effect_allele'))

    def test_empty_native_table_retains_normalized_schema(self):
        result=normalize_native(self.native_row().iloc[:0])
        self.assertEqual(len(result),0)
        self.assertTrue({'chrom_y','pos_y','ea_y','oa_y','beta_y','se_y','n_y'}.issubset(result.columns))

    def test_inverse_four_orientations(self):
        for code,ea,oa,b in [(10,'A','C',.2),(11,'C','A',-.2),(12,'T','G',.2),(13,'G','T',-.2)]:
            self.assertEqual(restore_original(ea,oa,b,code),('A','C',.2,'documented_inverse'))
        self.assertTrue(np.isnan(restore_original('A','T',.2,9)[2]))

    def test_palindrome_recovery_and_exact_variant_frequency(self):
        y=pd.DataFrame([dict(hm_code=6,effect_allele='T',other_allele='A',beta=-.3,
            chromosome='1',base_pair_location=210,variant_id='rs1',rsid='rs1',
            phenotype='p',standard_error=.1,neg_log10_p_value=3,p_value=.001,n_y=10,
            n_y_status='verified',source_y='official')])
        inst=pd.DataFrame([dict(chrom=1,pos38=210,pos37=100,variant='rs1')])
        af=pd.DataFrame([dict(chr=1,pos37=100,variant_y='rs1',ea_y='A',oa_y='T',eaf_y=.2,source_y='AF')])
        r=normalize_ebi(y,inst,af).iloc[0]
        self.assertEqual((r.ea_y,r.oa_y,r.beta_y,r.eaf_y),('A','T',.3,.2))
        self.assertEqual((r.build_y,r.pos_y),('GRCh38',210))
        y['variant_id']='rs2';y['rsid']='rs2'
        self.assertTrue(np.isnan(normalize_ebi(y,inst,af).iloc[0].eaf_y))


if __name__=='__main__':unittest.main()

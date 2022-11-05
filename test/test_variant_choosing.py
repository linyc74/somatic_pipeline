from somatic_pipeline.variant_choosing import VariantChoosing, BuildHeaderContigLines
from .setup import TestCase


class TestVariantChoosing(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = VariantChoosing(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa',
            vcfs=[
                f'{self.indir}/lofreq.vcf',
                f'{self.indir}/haplotype-caller.vcf',
                f'{self.indir}/vardict.vcf',
            ],
            min_num_callers=2,
        )
        expected = f'{self.indir}/variants.vcf'
        self.assertFileEqual(expected, actual)


class TestBuildHeaderContigLines(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = BuildHeaderContigLines(self.settings).main(
            ref_fa=f'{self.indir}/tiny.fa')
        expected = f'''\
##contig=<ID=chrA,length=10>
##contig=<ID=chrB,length=9>'''
        self.assertEqual(expected, actual)

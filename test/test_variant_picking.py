from somatic_pipeline.variant_picking import VariantPicking, BuildHeaderContigLines, GetChromToOrder
from .setup import TestCase


class TestVariantPicking(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = VariantPicking(self.settings).main(
            ref_fa=f'{self.indir}/chr9.fa',
            vcfs=[
                f'{self.indir}/lofreq.vcf',
                f'{self.indir}/muse.vcf',
                f'{self.indir}/mutect2.vcf',
                f'{self.indir}/somatic-sniper.vcf',
                f'{self.indir}/vardict.vcf',
                f'{self.indir}/varscan.vcf',
            ],
            min_snv_callers=3,
            min_indel_callers=1,
        )
        self.assertEqual(f'{self.workdir}/picked-variants.vcf', actual)
        self.assertFileEqual(f'{self.indir}/picked-variants.vcf', actual)


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


class TestGetChromToOrder(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = GetChromToOrder(self.settings).main(
            ref_fa=f'{self.indir}/tiny.fa')
        expected = {'chrA': 0, 'chrB': 1}
        self.assertEqual(expected, actual)

from somatic_pipeline.variant_filtering import FlagVariants, flag_variant, parse_criterion, RemoveVariants
from .setup import TestCase


class TestFlagVariants(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_mutect2(self):
        actual = FlagVariants(self.settings).main(
            vcf=f'{self.indir}/mutect2.vcf',
            variant_flagging_criteria='DP_25: 25<=DP<=25'
        )
        expected = f'{self.workdir}/mutect2-flagged.vcf'
        self.assertFileExists(expected, actual)

    def test_tiny(self):
        actual = FlagVariants(self.settings).main(
            vcf=f'{self.indir}/tiny.vcf',
            variant_flagging_criteria='LOW_DP: DP<20, HIGH_MQ: MQ>=30'
        )
        expected = f'{self.indir}/tiny-flagged.vcf'
        self.assertFileEqual(expected, actual)


class TestFunctions(TestCase):

    def test_flag_variant(self):
        variant = {'FILTER': '.', 'INFO': 'DP=1'}
        flag = 'FLAG'
        for criterion, expected_flag in [
            ('DP<2', flag),
            ('DP<=2', flag),
            ('1<DP<2', '.'),
            ('1<=DP<2', flag),
            ('1<DP<=2', '.'),
            ('1<=DP<=2', flag),
            ('DP>1', '.'),
            ('DP>=1', flag),
            ('2>DP>1', '.'),
            ('2>=DP>1', '.'),
            ('2>DP>=1', flag),
            ('2>=DP>=1', flag),
        ]:
            with self.subTest(criterion=criterion, expected_flag=expected_flag):
                flagged = flag_variant(
                    variant=variant,
                    flag=flag,
                    criterion=parse_criterion(criterion)
                )
                self.assertEqual(expected_flag, flagged['FILTER'])

    def test_add_second_flag(self):
        variant = {'FILTER': 'PASS', 'INFO': 'DP=1;MQ=15'}
        variant = flag_variant(
            variant=variant,
            flag='LOW_DP',
            criterion=parse_criterion('DP<10')
        )
        variant = flag_variant(
            variant=variant,
            flag='MID_MQ',
            criterion=parse_criterion('10<MQ<20')
        )
        expected = {'FILTER': 'PASS;LOW_DP;MID_MQ', 'INFO': 'DP=1;MQ=15'}
        self.assertEqual(expected, variant)

    def test_parse_criterion(self):
        for criterion, expected in [
            ('DP<2',     "Criterion(key='DP', range=(-inf, 2.0), equal_max=False, equal_min=False)"),
            ('DP<=2',    "Criterion(key='DP', range=(-inf, 2.0), equal_max=True, equal_min=False)"),
            ('1<DP<2',   "Criterion(key='DP', range=(1.0, 2.0), equal_max=False, equal_min=False)"),
            ('1<=DP<2',  "Criterion(key='DP', range=(1.0, 2.0), equal_max=False, equal_min=True)"),
            ('1<DP<=2',  "Criterion(key='DP', range=(1.0, 2.0), equal_max=True, equal_min=False)"),
            ('1<=DP<=2', "Criterion(key='DP', range=(1.0, 2.0), equal_max=True, equal_min=True)"),
            ('DP>1',     "Criterion(key='DP', range=(1.0, inf), equal_max=False, equal_min=False)"),
            ('DP>=1',    "Criterion(key='DP', range=(1.0, inf), equal_max=False, equal_min=True)"),
            ('2>DP>1',   "Criterion(key='DP', range=(1.0, 2.0), equal_max=False, equal_min=False)"),
            ('2>=DP>1',  "Criterion(key='DP', range=(1.0, 2.0), equal_max=True, equal_min=False)"),
            ('2>DP>=1',  "Criterion(key='DP', range=(1.0, 2.0), equal_max=False, equal_min=True)"),
            ('2>=DP>=1', "Criterion(key='DP', range=(1.0, 2.0), equal_max=True, equal_min=True)"),
        ]:
            with self.subTest(criterion=criterion, expected=expected):
                actual = parse_criterion(s=criterion)
                self.assertEqual(expected, repr(actual))


class TestRemoveVariants(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = RemoveVariants(self.settings).main(
            vcf=f'{self.indir}/tiny-flagged.vcf',
            flags=['LOW_DP']
        )
        expected = f'{self.indir}/tiny-flagged-variant-removal.vcf'
        self.assertFileEqual(expected, actual)

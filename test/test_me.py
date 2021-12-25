from .setup import TestCase
from gatk_pipeline.trimming import TrimGalore


class MyTest(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_me(self):
        TrimGalore(self.settings).main(
        )

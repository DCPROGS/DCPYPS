import os

from dcpyps.ekdist import ekscn

class TestSCNFileLoading:
    def setUp(self):
        self.infile = "../DCPYPS/dcpyps/tests/AChsim.scn"
        self.header = ekscn.read_header(self.infile)
        self.itint, self.iampl, self.iprops = ekscn.read_data(self.infile, self.header)
        
    def test_infile_exists(self):
        assert os.path.isfile(self.infile)
        
    def test_SCN_header_loading(self):
        assert self.header
        
    def test_interval_number(self):
        assert self.header['nint'] == 13948
        
    def test_interval_loading(self):    
        assert self.header['nint'] == len(self.itint)
        
    def test_flags(self):
        assert not self.iprops.any()
        
    def test_amplitudes(self):
        assert self.iampl[0] == 0.
        assert self.iampl[1] == 6.
        assert self.iampl[2] == 0.

    def tearDown(self):
        self.header = None
        self.itint, self.iampl, self.iprops = None, None, None
        


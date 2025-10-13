import pyIMSRG 

ms=pyIMSRG.ModelSpace(2,'He4','He4')
ut=pyIMSRG.UnitTest(ms)
passed=ut.SanityCheck()
try:
    assert passed
    print("All tests passed")
except AssertionError:  
    print("Tests failed")
    raise AssertionError

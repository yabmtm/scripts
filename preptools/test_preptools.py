### Testing script to be used with
### $ pytest -v

import os
import preptools

# define the test functions

def test_amber_installation():
    """Tests to see if AmberTools is currectly installed."""

    # does the standard install location for AMBERHOME exist?
    AMBERHOME = os.environ['AMBERHOME']
    print('Environment variable AMBERHOME = %s exists.'%AMBERHOME)

    # does tleap exist?
    tleap_path = os.path.join(AMBERHOME, 'bin/tleap')
    assert os.path.exists(tleap_path)

    # does antechamber exist?
    antechamber_path = os.path.join(AMBERHOME, 'bin/antechamber') 
    assert os.path.exists(antechamber_path)
 

# Run the tests!

test_amber_installation()

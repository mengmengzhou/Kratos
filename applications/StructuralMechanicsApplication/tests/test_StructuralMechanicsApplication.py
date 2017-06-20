# import Kratos
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication 
import KratosMultiphysics.SolidMechanicsApplication 
import KratosMultiphysics.StructuralMechanicsApplication 

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import SimpleMeshMovingTest as TSimpleMeshMovingTest
from SmallTests import DynamicBossakTests as TDynamicBossakTests
from SmallTests import DynamicNewmarkTests as TDynamicNewmarkTests
from SmallTests import SprismMembranePatchTests as TSprismMembranePatchTests
from SmallTests import SprismBendingPatchTests as TSprismBendingPatchTests
from SmallTests import ShellQ4ThinTensionTests as TShellQ4ThinTensionTests
from SmallTests import ShellQ4ThinBendingRollUpTests as TShellQ4ThinBendingRollUpTests
from SmallTests import ShellQ4ThinDrillingRollUpTests as TShellQ4ThinDrillingRollUpTests
from SmallTests import ShellQ4ThickBendingRollUpTests as TShellQ4ThickBendingRollUpTests
from SmallTests import ShellQ4ThickDrillingRollUpTests as TShellQ4ThickDrillingRollUpTests
from SmallTests import ShellT3ThinBendingRollUpTests as TShellT3ThinBendingRollUpTests
from SmallTests import ShellT3ThinDrillingRollUpTests as TShellT3ThinDrillingRollUpTests
from SmallTests import EigenQ4Thick2x2PlateTests as TEigenQ4Thick2x2PlateTests
from SmallTests import EigenTL3D8NCubeTests as TEigenTL3D8NCubeTests
from SmallTests import ShellThinQ4MembraneTest as TShellThinQ4MembraneTest
# ShellThickElement3D3N tests
from SmallTests import ShellThickElement3D3NLinearStaticTests as TShellT3ThickLinearStaticTests
from SmallTests import ShellThickElement3D3NNonLinearStaticTests as TShellT3ThickNonLinearStaticTests
from SmallTests import ShellThickElement3D3NLinearDynamicTests as TShellT3ThickLinearDynamicTests
from SmallTests import ShellThickElement3D3NNonLinearDynamicTests as TShellT3ThickNonLinearDynamicTests
# ShellThinElement3D4N tests
from SmallTests import ShellThinElement3D4NLinearStaticTests as TShellQ4ThinLinearStaticTests
from SmallTests import ShellThinElement3D4NNonLinearStaticTests as TShellQ4ThinNonLinearStaticTests
from SmallTests import ShellThinElement3D4NLinearDynamicTests as TShellQ4ThinLinearDynamicTests
from SmallTests import ShellThinElement3D4NNonLinearDynamicTests as TShellQ4ThinNonLinearDynamicTests

## NIGTHLY TESTS
# Shell test
from NightlyTests import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests
# CL tests
from NightlyTests import IsotropicDamageSimoJuPSTest    as TIsotropicDamageSimoJuPSTest

## VALIDATION TESTS
# SPRISM tests
from ValidationTests import SprismPanTests              as TSprismPanTests
# Eigenvalues tests
from ValidationTests import Eigen3D3NThinCircleTests    as TEigen3D3NThinCircleTests

def AssambleTestSuites():
    ''' Populates the test suites to run.
    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"
    Return
    ------
    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
#    # Basic moving mesh test
#    smallSuite.addTest(TSimpleMeshMovingTest('test_execution'))
#    # Dynamic basic tests
#    smallSuite.addTest(TDynamicBossakTests('test_execution'))
#    smallSuite.addTest(TDynamicNewmarkTests('test_execution'))
#    # Patch test Small Displacements
#    smallSuite.addTest(TSDTwoDShearQuaPatchTest('test_execution'))
#    smallSuite.addTest(TSDTwoDShearTriPatchTest('test_execution'))
#    smallSuite.addTest(TSDTwoDTensionQuaPatchTest('test_execution'))
#    smallSuite.addTest(TSDTwoDTensionTriPatchTest('test_execution'))
#    smallSuite.addTest(TSDThreeDShearHexaPatchTest('test_execution'))
#    smallSuite.addTest(TSDThreeDShearTetraPatchTest('test_execution'))
#    smallSuite.addTest(TSDThreeDTensionHexaPatchTest('test_execution'))
#    smallSuite.addTest(TSDThreeDTensionTetraPatchTest('test_execution'))
#    # Patch test Total Lagrangian
#    smallSuite.addTest(TTLTwoDShearQuaPatchTest('test_execution'))
#    smallSuite.addTest(TTLTwoDShearTriPatchTest('test_execution'))
#    smallSuite.addTest(TTLTwoDTensionQuaPatchTest('test_execution'))
#    smallSuite.addTest(TTLTwoDTensionTriPatchTest('test_execution'))
#    smallSuite.addTest(TTLThreeDShearHexaPatchTest('test_execution'))
#    smallSuite.addTest(TTLThreeDShearTetraPatchTest('test_execution'))
#    smallSuite.addTest(TTLThreeDTensionHexaPatchTest('test_execution'))
#    smallSuite.addTest(TTLThreeDTensionTetraPatchTest('test_execution'))
#    # Patch test Updated Lagrangian
#    smallSuite.addTest(TULTwoDShearQuaPatchTest('test_execution'))
#    smallSuite.addTest(TULTwoDShearTriPatchTest('test_execution'))
#    smallSuite.addTest(TULTwoDTensionQuaPatchTest('test_execution'))
#    smallSuite.addTest(TULTwoDTensionTriPatchTest('test_execution'))
#    smallSuite.addTest(TULThreeDShearHexaPatchTest('test_execution'))
#    smallSuite.addTest(TULThreeDShearTetraPatchTest('test_execution'))
#    smallSuite.addTest(TULThreeDTensionHexaPatchTest('test_execution'))
#    smallSuite.addTest(TULThreeDTensionTetraPatchTest('test_execution'))
#    # SPRISM tests
#    smallSuite.addTest(TSprismMembranePatchTests('test_execution'))
#    smallSuite.addTest(TSprismBendingPatchTests('test_execution'))
#    smallSuite.addTest(TShellQ4ThinTensionTests('test_execution'))
#    smallSuite.addTest(TShellQ4ThinBendingRollUpTests('test_execution'))
#    smallSuite.addTest(TShellQ4ThinDrillingRollUpTests('test_execution'))
#    smallSuite.addTest(TShellQ4ThickBendingRollUpTests('test_execution'))
#    smallSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
#    smallSuite.addTest(TShellT3ThinBendingRollUpTests('test_execution'))
#    smallSuite.addTest(TShellT3ThinDrillingRollUpTests('test_execution'))
#    # Eigenvalues tests
#    smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
#    smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))
#    smallSuite.addTest(TShellThinQ4MembraneTest('test_execution'))    
    # ShellThickElement3D3N tests
    smallSuite.addTest(TShellT3ThickLinearStaticTests('test_execution'))
    smallSuite.addTest(TShellT3ThickNonLinearStaticTests('test_execution'))
    smallSuite.addTest(TShellT3ThickLinearDynamicTests('test_execution'))
    smallSuite.addTest(TShellT3ThickNonLinearDynamicTests('test_execution'))
    # ShellThinElement3D4N tests
    smallSuite.addTest(TShellQ4ThinLinearStaticTests('test_execution'))
    smallSuite.addTest(TShellQ4ThinNonLinearStaticTests('test_execution'))
    smallSuite.addTest(TShellQ4ThinLinearDynamicTests('test_execution'))
    smallSuite.addTest(TShellQ4ThinNonLinearDynamicTests('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    # Shell tests
    nightSuite.addTest(TShellT3IsotropicScordelisTests('test_execution'))
    # CL tests
    nightSuite.addTest(TIsotropicDamageSimoJuPSTest('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    # SPRISM tests
    validationSuite.addTest(TSprismPanTests('test_execution'))
    # Eigenvalues tests
    validationSuite.addTest(TEigen3D3NThinCircleTests('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
#            TSimpleMeshMovingTest,
#            TDynamicBossakTests,
#            TDynamicNewmarkTests,
#            TSDTwoDShearQuaPatchTest,
#            TSDTwoDShearTriPatchTest,
#            TSDTwoDTensionQuaPatchTest,
#            TSDTwoDTensionTriPatchTest,
#            TSDThreeDShearHexaPatchTest,
#            TSDThreeDShearTetraPatchTest,
#            TSDThreeDTensionHexaPatchTest,
#            TSDThreeDTensionTetraPatchTest,
#            TTLTwoDShearQuaPatchTest,
#            TTLTwoDShearTriPatchTest,
#            TTLTwoDTensionQuaPatchTest,
#            TTLTwoDTensionTriPatchTest,
#            TTLThreeDShearHexaPatchTest,
#            TTLThreeDShearTetraPatchTest,
#            TTLThreeDTensionHexaPatchTest,
#            TTLThreeDTensionTetraPatchTest,
#            TULTwoDShearQuaPatchTest,
#            TULTwoDShearTriPatchTest,
#            TULTwoDTensionQuaPatchTest,
#            TULTwoDTensionTriPatchTest,
#            TULThreeDShearHexaPatchTest,
#            TULThreeDShearTetraPatchTest,
#            TULThreeDTensionHexaPatchTest,
#            TULThreeDTensionTetraPatchTest,
#            TSprismMembranePatchTests,
#            TSprismBendingPatchTests,
#            TShellQ4ThickBendingRollUpTests,
#            TShellQ4ThinTensionTests,
#            TShellQ4ThinBendingRollUpTests,
#            TShellQ4ThinDrillingRollUpTests,
#            TShellQ4ThickDrillingRollUpTests,
#            TShellT3ThinBendingRollUpTests,
#            TShellT3ThinDrillingRollUpTests,
#            TShellT3IsotropicScordelisTests,
#            TShellThinQ4MembraneTest,
            TShellT3ThickLinearStaticTests,
            TShellT3ThickNonLinearStaticTests,
            TShellT3ThickLinearDynamicTests,
            TShellT3ThickNonLinearDynamicTests,
            TShellQ4ThinLinearStaticTests,
            TShellQ4ThinNonLinearStaticTests,
            TShellQ4ThinLinearDynamicTests,
            TShellQ4ThinNonLinearDynamicTests
            ######TSprismPanTests
        ])
    )
    
    if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TEigenQ4Thick2x2PlateTests,
                TEigenTL3D8NCubeTests,
                TEigen3D3NThinCircleTests
            ])
        )
    else:
        print("FEASTSolver solver is not included in the compilation of the External Solvers Application")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
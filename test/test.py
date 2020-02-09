import sys
sys.path.insert(0, "..")

def test_generation():
    import FBMN_Evaluator

    test_filename = "networks/quantification_table-00001.csv"
    output_filename = "test_png.png"

    FBMN_Evaluator.fbmn_evaluate(test_filename, output_filename)

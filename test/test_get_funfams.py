# import os
# import pandas as pd
# import scripts.get_funfams as ff
#
#
# # Tests to check the validation of IDs
# def test_get_funfams():
#     df = pd.read_csv(
#         "test/files/funfams/baier_and_tokuriki_24_structures_24_generic_annotated.csv"
#     )
#
#     outpath = "test/files/funfams/"
#
#     if os.path.exists(outpath):
#         os.remove(outpath)
#
#     ff_df = ff.get_funfams(df, outpath)
#
#     print (ff_df['funfam_descriptions'])
#
#

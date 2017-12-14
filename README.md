# 使い方
* KEGG Pathwayを指定して含まれるタンパクの構造とinteractionを取得

'''
$ python main.py [inputfile] [out1] [out2]
'''

inputfileには `hsa04620` などと指定する．
out1には `pathway, kegg_id, uniprot_id, pdbid` が記述されて，out2には相互作用のある配列の対が記されます．

# analysisディレクトリ
PPI解析用のスクリプトがいくつかあります
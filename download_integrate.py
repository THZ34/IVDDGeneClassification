# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com

import os
import re
import anndata as ad
import pandas as pd
import gzip
import numpy as np
import requests
import seaborn as sns
from biokit.data import get_gse_sampleinfo, download_gse
from biokit.data._geo import get_gse_info
from biokit.preprocessing import genename_version_convert, detect_version, get_ensembl_genename_df, get_genename_df

import GEOparse
import mygene

# %%
gseids = ["GSE201396", "GSE131065", "GSE122429", "GSE186542", "GSE185728", "GSE113199", "GSE70362", "GSE242040",
          "GSE41883", "GSE27494", "GSE17077", "GSE147383"]

for gseid in gseids:
    os.makedirs(f'datasets/{gseid}', exist_ok=True)
    if not os.path.exists(f'datasets/{gseid}/{gseid}_sample_info.xlsx'):
        print(gseid)
        get_gse_sampleinfo(gseid).to_excel(f'datasets/{gseid}/{gseid}_sample_info.xlsx', index=False)
    # download_gse(gseid, output=f'datasets/{gseid}')

for gseid in gseids:
    if len(os.listdir(f'datasets/{gseid}')) == 1:
        print(gseid)

# %%
exp_df_dict = {}
sample_info_dict = {}
adata = ad.read_h5ad('T0020.leiden.celltype.h5ad')

# %% 基因名版本转换
# get_genename_df(list(range(76, 111)), '/data/project/tanghongzhen/data/database/annotation/ensembl', database='ensembl',
#                 outdir='/data/project/tanghongzhen/data/database/annotation')

gencode_genename_df = pd.read_csv('/data/project/tanghongzhen/data/database/annotation/gencode/genename.csv',
                                  index_col=0)
gencode_geneid_df = pd.read_csv('/data/project/tanghongzhen/data/database/annotation/gencode/geneid.csv', index_col=0)
gencode_geneid_df.fillna('NaN.1', inplace=True)
gencode_geneid_df = gencode_geneid_df.applymap(lambda x: x.split('.')[0])
gencode_geneid_df.replace({'NaN': np.nan}, inplace=True)

ensembl_genename_df = pd.read_csv('/data/project/tanghongzhen/data/database/annotation/ensembl/genename.csv',
                                  index_col=0)
ensembl_geneid_df = pd.read_csv('/data/project/tanghongzhen/data/database/annotation/ensembl/geneid.csv', index_col=0)

genename_df = pd.concat([gencode_genename_df, ensembl_genename_df], axis=1, join='outer')
genename_df.fillna(0, inplace=True)

geneid_df = pd.concat([gencode_geneid_df, ensembl_geneid_df], axis=1, join='outer')

entrez_id_df = pd.read_csv('/data/project/tanghongzhen/data/database/annotation/Homo_sapiens.gene_info.gz', index_col=0,
                           sep='\t')
entrez_to_symbol = entrez_id_df.set_index('GeneID')['Symbol'].to_dict()


# %% 基因长度
def extract_gene_id(attr, gene_id='gene_id'):
    match = re.search(f'{gene_id} "([^"]+)"', attr)
    return match.group(1) if match else None


def extract_gene_length(gtf_file):
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                         names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame",
                                "attribute"])
    gtf_df['gene_id'] = gtf_df['attribute'].map(lambda x: extract_gene_id(x, 'gene_id'))
    gtf_df['gene_name'] = gtf_df['attribute'].map(lambda x: extract_gene_id(x, 'gene_name'))
    gtf_df['gene_type'] = gtf_df['attribute'].map(lambda x: extract_gene_id(x, 'gene_type'))
    gtf_df['length'] = gtf_df['end'] - gtf_df['start'] + 1
    gene_length_dict = gtf_df[gtf_df['feature'] == 'transcript'].set_index('gene_name')['length'].to_dict()
    return gene_length_dict


def calculate_tpm(df, length_col="length"):
    # 提取 count 列（排除 length）
    count_cols = df.columns.difference([length_col])
    # Step 1: count / length
    rpk = df[count_cols].div(df[length_col], axis=0)
    # Step 2: 每个样本的 RPK 总和
    rpk_sum = rpk.sum()
    # Step 3: TPM = RPK / sum(RPK) * 1e6
    tpm_df = rpk.div(rpk_sum, axis=1) * 1e6
    # 设置 index 为 gene_name（保持原样）
    tpm_df.index = df.index
    tpm_df.index.name = "gene_name"

    return tpm_df


# %% GSE201396
# exp_df = pd.read_csv(f'datasets/GSE201396/GSE201396_FPKM.tsv', sep='\t', index_col=0)
# exp_df.set_index('gene_short_name', inplace=True)
# exp_df.drop('locus', axis=1, inplace=True)
# exp_df = exp_df.groupby(exp_df.index).mean()
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE201396/GSE201396_tpm.xlsx')
exp_tpm_df = pd.read_excel('datasets/GSE201396/GSE201396_tpm.xlsx', index_col=0)

sample_info = pd.read_excel(f'datasets/GSE201396/GSE201396_sample_info.xlsx')
sample_info['Title'] = sample_info['Title'].str.replace(' ', '-')

exp_df_dict['GSE201396'] = exp_tpm_df
sample_info_dict['GSE201396'] = sample_info

# %% GSE131065
# exp_df = pd.read_csv(f'datasets/GSE131065/GSE131065_invivoNPC.gene.expression.tables.FPKM.txt.gz', sep='\t',
#                      index_col=0)
# exp_df.set_index('gene_short_name', inplace=True)
# exp_df = exp_df.groupby(exp_df.index).mean()
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE131065/GSE131065_tpm.xlsx')
exp_tpm_df = pd.read_excel('datasets/GSE131065/GSE131065_tpm.xlsx', index_col=0)

sample_info = pd.read_excel(f'datasets/GSE131065/GSE131065_sample_info.xlsx')
sample_info['Title'].replace({'sample 124': 'CDS_124', 'sample135a': 'CDS_135_L1_L2', 'sample 135b': 'CDS_135_T12_L1',
                              'sample 147': 'CDS_147_L1_L2'}, inplace=True)

exp_df_dict['GSE131065'] = exp_tpm_df
sample_info_dict['GSE131065'] = sample_info

# %% GSE122429 (细胞实验,通过多能干细胞诱导产生NP细胞)
# exp_df = pd.read_csv(f'datasets/GSE122429/GSE122429_NPC.gene.expression.tables.FPKM.txt.gz', sep='\t', index_col=0)
# exp_df.set_index('gene_short_name', inplace=True)
# exp_df.drop(['tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'tss_id', 'locus'], axis=1, inplace=True)
# exp_df = exp_df.groupby(exp_df.index).mean()
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE122429/GSE122429_tpm.xlsx')
# exp_tpm_df = pd.read_excel('datasets/GSE122429/GSE122429_tpm.xlsx', index_col=0)
#
# sample_info = pd.read_excel(f'datasets/GSE122429/GSE122429_sample_info.xlsx')
# sample_info['Title'] = sample_info['Title'].map(lambda x: x.encode('latin1').decode('utf-8'))
# sample_info['Title'] = sample_info['Title'].str.replace('−', '-')  # U+2212数学负号换成连字符U+002D
# sample_info['Title'].replace({'ESC9': 'E9', 'ESC3': 'ES3', 'ImR90-iPSC': 'ImR90'}, inplace=True)
#
# exp_df_dict['GSE122429'] = exp_tpm_df
# sample_info_dict['GSE122429'] = sample_info

# %% GSE186542
# TODO: raw count, 需要先标准化到FPKM
# exp_df = pd.read_csv(f'datasets/GSE186542/GSE186542_gene_expression.txt.gz', sep='\t', index_col=0)
# exp_df.set_index('gene_name', inplace=True)
# exp_df.drop(['description', 'locus'], axis=1, inplace=True)
# exp_df = exp_df.groupby(exp_df.index).mean()
# exp_df.drop(
#     ['1-Dec', '1-Mar', '10-Mar', '11-Mar', '2-Mar', '3-Mar', '4-Mar', '5-Mar', '6-Mar', '7-Mar', '8-Mar', '9-Mar'],
#     axis=0, inplace=True)
# exp_df.columns = exp_df.columns.str.split('_count').str[0]
#
# current_version = detect_version(exp_df.index, genename_df)
#
# gene_length_dict = extract_gene_length(
#     '/data/project/tanghongzhen/data/database/annotation/ensembl/Homo_sapiens.GRCh38.94.gtf.gz')
# exp_df['length'] = exp_df.index.map(gene_length_dict)
# exp_df.dropna(subset='length', axis=0, inplace=True)
# exp_tpm_df = calculate_tpm(exp_df)
# exp_tpm_df.to_excel('datasets/GSE186542/GSE186542_tpm.xlsx')
exp_tpm_df = pd.read_excel('datasets/GSE186542/GSE186542_tpm.xlsx', index_col=0)

sample_info = pd.read_excel(f'datasets/GSE186542/GSE186542_sample_info.xlsx')

exp_df_dict['GSE186542'] = exp_tpm_df
sample_info_dict['GSE186542'] = sample_info

# %% GSE185728 (人NP细胞接种到不同刚度的水凝胶之后的转录组和表观组, 不是组织)
# exp_df = pd.read_csv(f'datasets/GSE185728/GSE185728_20210107_ZhouTaiFeng_mRNA-seq_hisat.counts.txt.gz', sep='\t',
#                      index_col=0, header=1)
# sample_info = pd.read_excel(f'datasets/GSE185728/GSE185728_sample_info.xlsx')

# %% GSE113199
# 样本名带CoPP的使用了CoPP药物处理，去掉
# TODO: raw count, 需要先标准化到FPKM
# exp_df = pd.read_csv(f'datasets/GSE113199/GSE113199_All_Counts.txt.gz', sep='\t', index_col=0)
# exp_df.drop(['GeneType'], axis=1, inplace=True)
# exp_df = exp_df.loc[:, ~exp_df.columns.str.contains('CoPP')]
#
# current_version = detect_version(exp_df.index, genename_df)[0]
# gene_length_dict = extract_gene_length(
#     '/data/project/tanghongzhen/data/database/annotation/ensembl/Homo_sapiens.GRCh38.85.gtf.gz')
# exp_df['length'] = exp_df.index.map(gene_length_dict)
# exp_df.dropna(subset='length', axis=0, inplace=True)
# exp_tpm_df = calculate_tpm(exp_df)
# exp_tpm_df.to_excel('datasets/GSE113199/GSE113199_tpm.xlsx')
exp_tpm_df = pd.read_excel('datasets/GSE113199/GSE113199_tpm.xlsx', index_col=0)

sample_info = pd.read_excel(f'datasets/GSE113199/GSE113199_sample_info.xlsx')
sample_info['Title'] = sample_info['Title'].map(lambda x: x.encode('latin1').decode('utf-8'))
sample_info = sample_info[sample_info['Title'].isin(exp_tpm_df.columns.to_list())]

exp_df_dict['GSE113199'] = exp_tpm_df
sample_info_dict['GSE113199'] = sample_info

# %% GSE70362
# TODO: 原始数据是CEL格式, 没有表达矩阵，需要先用affy处理成表达矩阵
# TODO: 其中包含纤维环和髓核样本, 需要去掉纤维环
sample_info = pd.read_excel('datasets/GSE70362/GSE70362_sample_info.xlsx')

# exp_df = pd.read_csv('datasets/GSE70362/GSE70362_raw_exp_df.csv', index_col=0)
# # 转换成基因名
# info_dict = {}
# get_gse_info('GSE70362', info_dict)
# gpl = GEOparse.get_GEO("GPL570", destdir="/data/project/tanghongzhen/data/.cache")
# gpl_table = gpl.table
# gpl_table['Gene Symbol'] = gpl_table['Gene Symbol'].str.split(' /// ').str[0]
# exp_df.index = exp_df.index.map(gpl_table.set_index('ID')['Gene Symbol'].to_dict())
# exp_df = exp_df[~exp_df.index.isna()]
# exp_df = exp_df.groupby(exp_df.index).mean()
#
# # gsmid转换样本名
# exp_df.columns = exp_df.columns.str.split('_').str[0]
# exp_df.columns = exp_df.columns.map(sample_info.set_index('gsmid')['Title'].to_dict())
#
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE70362/GSE70362_tpm.xlsx')
exp_tpm_df = pd.read_excel('datasets/GSE70362/GSE70362_tpm.xlsx', index_col=0)

# 去掉纤维环
sample_info = sample_info[sample_info['tissue'] != 'Annulus fibrosus']
exp_tpm_df = exp_tpm_df[sample_info['Title']]

exp_df_dict['GSE70362'] = exp_tpm_df
sample_info_dict['GSE70362'] = sample_info

# %% GSE242040 基因数太少
# # TODO: 基因名使用的是geneid, 需要转换成genename
exp_df = pd.read_csv(f'datasets/GSE242040/GSE242040_ExpressionTable.txt.gz', sep='\t', index_col=0)
# exp_df.dropna(inplace=True, axis=1)
#
# union_df = []
# for col in geneid_df.columns:
#     union_df.append([col, len(set(exp_df.index) & set(geneid_df[col]))])
# union_df = pd.DataFrame(union_df, columns=['version', 'count'])
#
# exp_df.index = exp_df.index.map()

# sample_info = pd.read_excel(f'datasets/GSE242040/GSE242040_sample_info.xlsx')

# %% GSE41883
# TODO: 原始数据是CEL格式, 没有表达矩阵，需要先用affy处理成表达矩阵
# 纤维环组织
# sample_info = pd.read_excel('datasets/GSE41883/GSE41883_sample_info.xlsx')

# exp_df = pd.read_csv('datasets/GSE41883/GSE41883_raw_exp_df.csv', index_col=0)
# # 转换成基因名
# gpl = GEOparse.get_GEO("GPL1352", destdir="/data/project/tanghongzhen/data/.cache")
# gpl_table = gpl.table
# gpl_table['Gene Symbol'] = gpl_table['Gene Symbol'].str.split(' /// ').str[0]
# exp_df.index = exp_df.index.map(gpl_table.set_index('ID')['Gene Symbol'].to_dict())
# exp_df = exp_df[~exp_df.index.isna()]
# exp_df = exp_df.groupby(exp_df.index).mean()
#
# # gsmid转换样本名
# exp_df.columns = exp_df.columns.str.split('_').str[0]
# exp_df.columns = exp_df.columns.map(sample_info.set_index('gsmid')['Title'].to_dict())
#
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE41883/GSE41883_tpm.xlsx')
# exp_tpm_df = pd.read_excel('datasets/GSE41883/GSE41883_tpm.xlsx', index_col=0)
#
# exp_df_dict['GSE41883'] = exp_tpm_df
# sample_info_dict['GSE41883'] = sample_info

# %% GSE27494
# TODO: 原始数据是CEL格式, 没有表达矩阵，需要先用affy处理成表达矩阵
# 纤维环组织
# sample_info = pd.read_excel('datasets/GSE27494/GSE27494_sample_info.xlsx')

# exp_df = pd.read_csv('datasets/GSE27494/GSE27494_raw_exp_df.csv', index_col=0)
# # 转换成基因名
# gpl = GEOparse.get_GEO("GPL1352", destdir="/data/project/tanghongzhen/data/.cache")
# gpl_table = gpl.table
# gpl_table['Gene Symbol'] = gpl_table['Gene Symbol'].str.split(' /// ').str[0]
# exp_df.index = exp_df.index.map(gpl_table.set_index('ID')['Gene Symbol'].to_dict())
# exp_df = exp_df[~exp_df.index.isna()]
# exp_df = exp_df.groupby(exp_df.index).mean()
#
# # gsmid转换样本名
# exp_df.columns = exp_df.columns.str.split('.').str[0]
# exp_df.columns = exp_df.columns.map(sample_info.set_index('gsmid')['Title'].to_dict())
#
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE27494/GSE27494_tpm.xlsx')
# exp_tpm_df = pd.read_excel('datasets/GSE27494/GSE27494_tpm.xlsx', index_col=0)
#
# exp_df_dict['GSE27494'] = exp_tpm_df
# sample_info_dict['GSE27494'] = sample_info

# %% GSE17077
# TODO: 原始数据是CEL格式, 没有表达矩阵，需要先用affy处理成表达矩阵
sample_info = pd.read_excel('datasets/GSE17077/GSE17077_sample_info.xlsx')

# exp_df = pd.read_csv('datasets/GSE17077/GSE17077_raw_exp_df.csv', index_col=0)
# # 转换成基因名
# gpl = GEOparse.get_GEO("GPL1352", destdir="/data/project/tanghongzhen/data/.cache")
# gpl_table = gpl.table
# gpl_table['Gene Symbol'] = gpl_table['Gene Symbol'].str.split(' /// ').str[0]
# exp_df.index = exp_df.index.map(gpl_table.set_index('ID')['Gene Symbol'].to_dict())
# exp_df = exp_df[~exp_df.index.isna()]
# exp_df = exp_df.groupby(exp_df.index).mean()
#
# # gsmid转换样本名
# exp_df.columns = exp_df.columns.str.split('.').str[0]
# exp_df.columns = exp_df.columns.map(sample_info.set_index('gsmid')['Title'].to_dict())
#
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE17077/GSE17077_tpm.xlsx')
exp_tpm_df = pd.read_excel('datasets/GSE17077/GSE17077_tpm.xlsx', index_col=0)

exp_df_dict['GSE17077'] = exp_tpm_df
sample_info_dict['GSE17077'] = sample_info

# %% GSE147383
# TODO: 原始数据是CEL格式, 没有表达矩阵，需要先用affy处理成表达矩阵
# TODO:有纤维环和髓核,去掉纤维环
sample_info = pd.read_excel('datasets/GSE147383/GSE147383_sample_info.xlsx')

# exp_df = pd.read_csv('datasets/GSE147383/GSE147383_raw_exp_df.csv', index_col=0)
# # # 转换成基因名
# gpl = GEOparse.get_GEO("GPL570", destdir="/data/project/tanghongzhen/data/.cache")
# gpl_table = gpl.table
# gpl_table['Gene Symbol'] = gpl_table['Gene Symbol'].str.split(' /// ').str[0]
# exp_df.index = exp_df.index.map(gpl_table.set_index('ID')['Gene Symbol'].to_dict())
# exp_df = exp_df[~exp_df.index.isna()]
# exp_df = exp_df.groupby(exp_df.index).mean()
#
# # gsmid转换样本名
# exp_df.columns = exp_df.columns.str.split('_').str[0]
# exp_df.columns = exp_df.columns.map(sample_info.set_index('gsmid')['Title'].to_dict())
#
# exp_tpm_df = exp_df / exp_df.sum() * 1e6
# exp_tpm_df.to_excel('datasets/GSE147383/GSE147383_tpm.xlsx')
exp_tpm_df = pd.read_excel('datasets/GSE147383/GSE147383_tpm.xlsx', index_col=0)

sample_info = sample_info[sample_info['tissue'] != 'annulus fibrosus']
exp_tpm_df = exp_tpm_df[sample_info['Title']]

exp_df_dict['GSE147383'] = exp_tpm_df
sample_info_dict['GSE147383'] = sample_info

# %%
exp_tpm_df = exp_df_dict['GSE201396']
genes = set(exp_tpm_df.index)
print(len(genes))
for gseid in ['GSE201396', "GSE131065", "GSE122429", "GSE186542", "GSE185728", "GSE113199", "GSE70362", "GSE242040",
              "GSE41883", "GSE27494", "GSE17077", "GSE147383"]:
    if gseid in exp_df_dict:
        exp_tpm_df = exp_df_dict[gseid]
        sample_info = sample_info_dict[gseid]
        sample_info['gseid'] = gseid
        # genes = genes & set(exp_tpm_df.index)
        # print(gseid, len(genes), exp_tpm_df.shape[0])
        # print(exp_tpm_df.index.is_unique)
        print(gseid, exp_tpm_df.columns.intersection(sample_info['Title']).shape[0],
              exp_tpm_df.columns.intersection(sample_info['Title']).shape[0] == exp_tpm_df.shape[1])

all_sample_exp_tpm_df = pd.concat([exp_df_dict[i] for i in exp_df_dict], axis=1, join='inner')
all_sample_info_df = pd.concat([sample_info_dict[i] for i in sample_info_dict], axis=0, join='outer')
all_sample_info_df.to_excel('datasets/all_sample_info.xlsx')
all_sample_exp_tpm_df.to_excel('datasets/all_sample_exp_tpm.xlsx')

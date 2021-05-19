#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

import pandas as pd
from loguru import logger

df = pd.read_csv('haploid.cell.line.expression.csv')
df = df.rename(columns={'gene': 'Gene', 'chrom': 'Chrom', 'start': 'Start', 'end': 'End'})
logger.info(df.shape)


def gene_expression():
    csv = 'gene.expression.csv'
    if os.path.isfile(csv):
        logger.info(f'File {csv} already exists, loading data ...')
        dd = pd.read_csv(csv)
    else:
        dd = pd.DataFrame()
        for cell_line in ('A11.C636', 'A11', 'E9', 'HAP1.C597', 'KBM7'):
            logger.info(f'Processing cell line {cell_line}')
            dx = df[['Gene', 'Transcript_ID'] + [cell_line]]
            dx = dx.sort_values(by=['Gene', 'Transcript_ID', cell_line], ascending=[True, True, False])
            dx = dx.drop_duplicates(subset=['Gene'])
            dx = dx.drop(columns=['Transcript_ID'])
            logger.info(f'dx: {dx.shape}')
            dd = dx if dd.empty else pd.merge(dd, dx, on=['Gene'])
            logger.info(f'dd: {dd.shape}')
            logger.info(f'dd\n{dd.head(3)}')
        logger.info(dd.shape)
        logger.info(dd.head())
        logger.info(f'Saving results ...')
        dd = dd.drop_duplicates()
        logger.info(dd.shape)
        dd.to_csv(csv, index=False, float_format='%.2f')
    return dd


def transcript_expression():
    csv = 'transcript.expression.csv'
    if os.path.isfile(csv):
        logger.info(f'File {csv} already exists, loading data ...')
        dd = pd.read_csv(csv)
    else:
        dd = df.drop(columns='Chrom,Start,End,strand,utr5_start,utr5_end,utr3_start,utr3_end'.split(','))
        dd.to_csv(csv, index=False, float_format='%.2f')
    return dd
    
    
def transcript_coordinates():
    csv = 'transcript.coordinates.csv'
    if os.path.isfile(csv):
        logger.info(f'File {csv} already exists, loading data ...')
        dc = pd.read_csv(csv)
    else:
        dc = df.drop(columns='A11.C636,A11,E9,HAP1.C597,KBM7'.split(','))
        dc = dc.drop_duplicates(subset=['Transcript_ID'])
        dc.to_csv(csv, index=False, float_format='%.2f')
    return dc
    
    
def intersect_rbp():
    rbp_gene_expression, rbp_transcript_expression = 'RBP_gene_expression.csv', 'RBP_transcript_expression.csv'
    rbp = pd.read_csv('RBP.list.tsv', sep='\t')
    columns = {'gene name': 'Gene', 'description': 'Description', 'gene id': 'ENSG',
               'protein id': 'ENSP', 'number of domains': 'Number_of_domains',
               'domains[count]': 'Domains_count', 'paralogous family ID': 'Paralogous_family_ID',
               'members of paralogous family ': 'Members_of_paralogus_family',
               'representative family member': 'Representative family member',
               'paralog [ % identity representative member:paralog , % identity paralog: representative member]':
                   'paralog',
               'consensus RNA target': 'Consensus_RNA_target',
               'putative RNA target': 'Putative_RNA_target',
               'supporting evidence (# pubmed ID)': 'PumMed_ID',
               'category in Baltz et al.': 'Category_in_Baltz_et_al.',
               'category in Castello et al.': 'Category_in_Castello_et_al.',
               'found in both human mass spec': 'Found_in_both_human_mass_spec',
               'human cell lines mass spec': 'Human_cell_lines_mass_spec',
               'category in Kwon et al. mouse mESC': 'category_in_Kwon_et_al._mouse_mESC',
               'Toronto RBPDB': 'Toronto_RBPDB',
               'GO RNA related category': 'GO_RNA_related_category',
               'GO RNA related category description': 'GO_RNA_related_category_description',
               'GO general category': 'GO_general_category',
               'GO general category description': 'GO general category description',
               'OMIM': 'OMIM'}
    rbp = rbp.rename(columns=columns)
    rbp[['Description', 'HGNC']] = rbp['Description'].str.split('\ \[Source:HGNC Symbol;Acc:', expand=True)
    rbp['HGNC'] = rbp['HGNC'].str.rstrip(']')
    print(rbp.head())
    print(list(rbp.columns))
    ge, te = gene_expression(), transcript_expression()
    if os.path.isfile(rbp_gene_expression):
        logger.info(f'File {rbp_gene_expression} already exists, loading data ...')
        dg = pd.read_csv(rbp_gene_expression)
    else:
        dg = pd.merge(rbp[['Gene', 'Description', 'HGNC']], ge, on='Gene')
        dg.to_csv(rbp_gene_expression, index=False)
    if os.path.isfile(rbp_transcript_expression):
        logger.info(f'File {rbp_transcript_expression} already exists, loading data ...')
        dt = pd.read_csv(rbp_transcript_expression)
    else:
        dt = pd.merge(rbp, te, on='Gene')
        dt.to_csv(rbp_transcript_expression, index=False)
    return dg, dt
        

if __name__ == '__main__':
    # gene_expression()
    # transcript_expression()
    # transcript_coordinates()
    intersect_rbp()

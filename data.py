#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from loguru import logger

df = pd.read_csv('haploid.cell.line.expression.csv')
df = df.rename(columns={'gene': 'Gene', 'chrom': 'Chrom', 'start': 'Start', 'end': 'End'})
logger.info(df.shape)


def gene_expression():
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
    dd.to_csv('gene.expression.csv', index=False, float_format='%.2f')


def transcript_expression():
    dd = df.drop(columns='Chrom,Start,End,strand,utr5_start,utr5_end,utr3_start,utr3_end'.split(','))
    dd.to_csv('transcript.expression.csv', index=False, float_format='%.2f')
    
    
def transcript_coordinates():
    dc = df.drop(columns='A11.C636,A11,E9,HAP1.C597,KBM7'.split(','))
    dc = dc.drop_duplicates(subset=['Transcript_ID'])
    dc.to_csv('transcript.coordinates.csv', index=False, float_format='%.2f')
    

if __name__ == '__main__':
    # gene_expression()
    # transcript_expression()
    transcript_coordinates()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
from bottle import route, run, template
import requests
from loguru import logger

logger.remove()
logger.add(sys.stderr, level="DEBUG", format="<light-green>[{time:HH:mm:ss}]</light-green> <level>{message}</level>")

cell_lines = ['A11.C636', 'A11', 'E9', 'HAP1.C597', 'KBM7']
rbp_columns = {'gene name': 'Gene', 'description': 'Description', 'gene id': 'ENSG',
               'protein id': 'ENSP', 'number of domains': 'Number of domains',
               'domains[count]': 'Domains count', 'paralogous family ID': 'Paralogous family ID',
               'members of paralogous family ': 'Members of paralogus family',
               'representative family member': 'Representative family member',
               'paralog [ % identity representative member:paralog , % identity paralog: representative member]': 'Paralog',
               'consensus RNA target': 'Consensus RNA target',
               'putative RNA target': 'Putative RNA target',
               'supporting evidence (# pubmed ID)': 'PubMed ID',
               'category in Baltz et al.': 'Category in Baltz et al.',
               'category in Castello et al.': 'Category in Castello et al.',
               'found in both human mass spec': 'Found in both human mass spec',
               'human cell lines mass spec': 'Human cell lines mass spec',
               'category in Kwon et al. mouse mESC': 'Category in Kwon et al. mouse mESC',
               'Toronto RBPDB': 'Toronto RBPDB',
               'GO RNA related category': 'GO RNA related category',
               'GO RNA related category description': 'GO RNA related category description',
               'GO general category': 'GO general category',
               'GO general category description': 'GO general category description',
               'OMIM': 'OMIM'}


def rbps():
    rbp = pd.read_csv('RBP.list.tsv', sep='\t')
    rbp = rbp.rename(columns=rbp_columns)
    rbp[['Description', 'HGNC']] = rbp['Description'].str.split('\ \[Source:HGNC Symbol;Acc:', expand=True)
    rbp['Description'] = rbp['Description'].str.replace('\ \[Source:.*', '', regex=True)
    rbp['HGNC'] = rbp['HGNC'].str.rstrip(']')
    return rbp


def expressions():
    df = pd.read_csv('haploid.cell.line.expression.csv')
    df = df.rename(columns={'gene': 'Gene', 'chrom': 'Chrom', 'start': 'Start', 'end': 'End'})
    return df


def rbp_expressions():
    rbp, dd = rbps(), expressions()
    df = pd.merge(rbp, dd, on='Gene')
    return df


def rbp_gene_expressions():
    csv = 'haploid.rbp.gene.expressions.csv'
    if os.path.isfile(csv):
        dd = pd.read_csv(csv)
        logger.info(f'Haploid RBP expressions already exist, loading data ...')
    else:
        df, dd = rbp_expressions(), pd.DataFrame()
        for cell_line in cell_lines:
            logger.info(f'Processing cell line {cell_line}')
            dx = df[['Gene', 'Transcript ID', 'Description', cell_line]]
            dx = dx.sort_values(by=['Gene', 'Transcript ID', cell_line], ascending=[True, True, False])
            dx = dx.drop_duplicates(subset=['Gene'])
            dx = dx.drop(columns=['Transcript ID'])
            dd = dx if dd.empty else pd.merge(dd, dx, on=['Gene', 'Description'], how='outer')
        logger.info(f'Saving results ...')
        dd = dd.drop_duplicates()
        dd = dd.fillna(0.0)
        dd.to_csv(csv, index=False, float_format='%.2f')
    return dd


def rbp_transcript_expressions():
    csv = 'haploid.rbp.transcript.expression.csv'
    if os.path.isfile(csv):
        logger.info(f'File {csv} already exists, loading data ...')
        dd = pd.read_csv(csv)
    else:
        rbp, df = rbps(), expressions()
        dd = pd.merge(rbp, df, on=['Gene'])
        dd.to_csv(csv, index=False, float_format='%.2f')
    return dd


def get_transcript(gene='', transcript=''):
    df = rbp_transcript_expressions()
    gt = df[df['Gene'] == gene].copy() if gene else df[df['Transcript ID'] == transcript].copy()
    gene, transcript = gt.iloc[0][['Gene', 'Transcript ID']]
    table = gt[['Transcript ID'] + cell_lines].copy()
    table['Transcript ID'] = [f'<a href="/transcript/{t}">{t}</a>' for t in table['Transcript ID']]
    table = table.to_html(index=False, escape=False, classes="table table-striped table-bordered",
                          table_id="transcript_table", justify='center')

    summary = gt.drop_duplicates(subset=['Gene']).copy()
    summary['Coordinate'] = [ucsc_url(row.Chrom, row.Start, row.End) for row in summary.itertuples()]
    columns = ['Gene', 'Description', 'Coordinate', 'Domains count', 'Paralogous family ID',
               'Members of paralogus family', 'Representative family member', 'Paralog', 'Consensus RNA target',
               'Putative RNA target', 'PubMed ID', 'GO RNA related category', 'GO RNA related category description',
               'GO general category', 'GO general category description', 'OMIM']
    summary_title = f'{gene} ({summary.iloc[0]["ENSG"]})'
    summary = summary[columns]
    summary = summary.set_index('Gene').transpose()
    summary['Gene'] = summary.index
    summary = summary[['Gene', gene]]
    summary = summary.to_html(escape=False, header=False, classes="table table-striped table-bordered",
                              table_id="gene_description", index=False)
    summary = summary.replace('<tbody>', """  <colgroup>
        <col style="width:30%; text-align:right">
        <col style="width:70%; text-align:left">
        </colgroup>
        <tbody>""")
    return gene, transcript, gt, table, summary_title, summary


def get_sequence(chrom, start, end, strand, utr5start, utr5end, utr3start, utr3end):
    server = "http://grch37.rest.ensembl.org"
    # start, end = (utr3start, utr3end) if strand == '+' else (utr3start, utr5end)
    ext = f"/sequence/region/human/{chrom.replace('chr', '')}:{start}..{end}:{1 if strand == '+' else -1}?"
    r = requests.get(server + ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    seq = r.text
    length = len(seq)
    utr5_length, utr3_length = utr5end - utr5start + 1, utr3end - utr3start + 1
    if strand == '+':
        start_codon_start_pos, stop_codon_start_pos = utr5end + 2, utr3start - 2
        start_codon_start, stop_codon_start = utr5_length + 1, length - utr3_length - 2
    else:
        start_codon_start_pos, stop_codon_start_pos = utr5start, utr3end + 4
        start_codon_start, stop_codon_start = utr5_length - 4, length - utr3_length + 3
    sequence = [
        ('utr5', seq[:start_codon_start]),
        ('start_codon', seq[start_codon_start: start_codon_start + 3]),
        ('cds', seq[start_codon_start + 3: stop_codon_start]),
        ('stop_codon', seq[stop_codon_start: stop_codon_start + 3]),
        ('utr3', seq[stop_codon_start + 3:])
    ]
    a = (f'<a class="btn btn-outline-success btn-sm" '
         f'href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={chrom}%3A{start}%2D{end}" '
         f'target="_blank">{chrom}:{start:,}-{end:,}<i class="fas fa-external-link-alt pl-2"></i></a>')
    sequence = '\n'.join([f'<span class="{c}">{t}</span>' for c, t in sequence])
    return sequence, a


@route('/')
def home(number=0):
    df = rbp_gene_expressions()
    df['Gene'] = [f'<a href="/gene/{gene}">{gene}</a>' for gene in df['Gene']]
    table = df.to_html(index=False, classes="table table-striped table-bordered", table_id="gene_expression",
                       justify='center', escape=False)
    return template('index.html', table=table)


@route('/gene/<gene>')
def gene_transcript(gene):
    _, _, _, table, summary_title, summary = get_transcript(gene=gene)
    return template('transcript.html', table=table, summary_title=summary_title, summary=summary,
                    table_title='', sequence_header='', sequence='')


@route('/transcript/<transcript>')
def transcript_expression(transcript):
    gene, _, gt, table, summary_title, summary = get_transcript(transcript=transcript)
    t = gt[gt['Transcript ID'] == transcript].iloc[0]
    sequence, a = get_sequence(t.Chrom, t.Start, t.End, t.strand, t.utr5_start, t.utr5_end, t.utr3_start, t.utr3_end)
    sequence_header = f'> {transcript} ({gene}|{t.ENSG}) {a}'
    return template('transcript.html', table=table, summary_title=summary_title, summary=summary,
                    sequence_header=sequence_header, sequence=sequence)


def ucsc_url(chrom, start, end, strand=''):
    s = f'{chrom}:{start}-{end}'
    if strand:
        s = f'{s}:{strand}'
    url = (f'<a href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={chrom}%3A{start}%2D{end}" '
           f'target="_blank">{s}</a>')
    return url


if __name__ == '__main__':
    run(host='localhost', port=8080, debug=True)
